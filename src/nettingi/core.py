
import numpy as np
import os
import psutil

import math as math

import time

from blimpy import GuppiRaw

from .utils import (
    template_check_outfile,
    template_check_nblocks,
    template_print_header,
    template_save_npy,
    template_calc_ave,
    template_print_flagstats,
    repl_nans_jit,
    repl_zeros,
    statistical_noise_fir,
    template_guppi_format,
)
from .reduction import raw2spec_god
from .plotting import load_raw_flags


class mitigateRFI:
    """
    blah blah blah
    
    Parameters
    ----------
        det_method : str
            RFI detection method used.
        inputfile : str
            input file
        repl_method : str
            replacement method
        cust : str
            
    
    """
    def __init__(self):


        #default/hardcoded attributes
        self.in_dir = '/jetstor/scratch/rfimit/unmitigated/rawdata/'


    def run_all(self):
        #do all the rfi mitigation steps


        pp = psutil.Process(os.getpid())

        start_time = time.time()
        if self.output_bool:
            template_check_outfile(self.infile_raw_full,self.outfile_raw_full)
            out_rawFile = open(self.outfile_raw_full,'rb+')

        
        template_check_nblocks(self._rawFile,self.mb)
        numblocks = self._rawFile.find_n_data_blocks()

        for bi in range(numblocks//self.mb):
            print('------------------------------------------')
            print(f'Block: {(bi*self.mb)+1}/{numblocks}')

            bstart = time.time()

            #print header for the first block
            if bi == 0:
                headersize = template_print_header(self._rawFile)


            #loading multiple blocks at once?        
            for mb_i in range(self.mb):
                if mb_i==0:
                    header,data = self._rawFile.read_next_data_block()
                    data = np.copy(data)
                    d1s = data.shape[1]
                else:
                    h2,d2 = self._rawFile.read_next_data_block()
                    data = np.append(data,np.copy(d2),axis=1)
            #data = np.ascontiguousarray(data)

            #find data shape
            # print(f'Data shape: {data.shape} || block size: {data.nbytes}')

            #save raw data?
            if self.rawdata:
                template_save_npy(data,bi,self.npybase)    
       

            spect_block = template_calc_ave(data, self.ave_factor)


            #===============================================
            #***********************************************

            if self.det_method == 'SK':
                print('SK mitigation')
                flags_block, ss_sk_block, ms_sk_block = self.SK_detection(data)
                if bi == 0:
                    self.ss_sk_all = ss_sk_block
                    self.ms_sk_all = ms_sk_block
                else:
                    self.ss_sk_all = np.concatenate((self.ss_sk_all, ss_sk_block),axis=1)
                    self.ms_sk_all = np.concatenate((self.ms_sk_all, ms_sk_block),axis=1)

            elif self.det_method == 'IQRM':
                flags_block = self.iqrm_detection(data)

            elif self.det_method == 'AOF':
                flags_block = self.aof_detection(data)

            elif self.det_method == 'MAD':
                flags_block = self.mad_detection(data)
        
            elif self.det_method == 'SE':
                flags_block, zsc_block = self.se_detection(data)
                if bi == 0:
                    self.zsc_all = zsc_block
                else:
                    self.zsc_all = np.concatenate((self.zsc_all, zsc_block),axis=1)

            elif self.det_method == 'Conv':
                flags_block = self.conv_detection(data)

            #***********************************************
            #===============================================



            #track intermediate results
            if bi == 0:
                self.flags_all = flags_block
                self.spect_all = spect_block
            else:
                self.flags_all = np.concatenate((self.flags_all, flags_block),axis=1)
                self.spect_all = np.concatenate((self.spect_all, spect_block),axis=1)

            if (self.det_method == 'AOF') or (self.det_method == 'MAD'):
                block_fname = str(bi).zfill(3)
                save_fname = self.output_mit_srdp_dir+self.npybase+'_flags_block'+block_fname+'.npy'
                np.save(save_fname,flags_block)
                self.flags_all = np.empty((data.shape[0],1,data.shape[2]))




             #print(f'MEM: spect: {self.spect_all.nbytes/1e9} // flags: {self.flags_all.nbytes/1e9}')
            mu = pp.memory_info()
            print(f'Total RAM usage: {mu[0]/2.**30} GB')
            #track flags


            template_print_flagstats(flags_block, False)


            #now flag shape is (chan,spectra,pol)
            #apply union of flags between the pols
            if True:
                flags_block[:,:,0][flags_block[:,:,1]==1]=1
                flags_block[:,:,1][flags_block[:,:,0]==1]=1


            ts_factor = data.shape[1] // flags_block.shape[1]
            rlen = flags_block.shape[1] * ts_factor
                

            if self.repl_method == 'nans':
                data[:,:rlen,:] = repl_nans_jit(data[:,:rlen,:],flags_block)

            if self.repl_method == 'zeros':
                #replace data with zeros
                data[:,:rlen,:] = repl_zeros(data[:,:rlen,:],flags_block)

            if self.repl_method == 'previousgood':
                #replace data with previous (or next) good
                print('Cannot replace with previous good data currently, switching to statistical noise')
                data[:,:rlen,:] = statistical_noise_fir(data[:,:rlen,:],flags_block,ts_factor)

            if self.repl_method == 'stats':
                #replace data with statistical noise derived from good datapoints
                data[:,:rlen,:] = statistical_noise_fir(data[:,:rlen,:],flags_block,ts_factor)


            #save the regenerated spectra
            regen_block = template_calc_ave(data, self.ave_factor)

            if bi == 0:
                self.regen_all = regen_block
            else:
                self.regen_all = np.concatenate((self.regen_all, regen_block),axis=1)


            #write back raw data
            if self.output_bool:

                #print('Re-formatting data and writing back to file...')

                for mb_i in range(self.mb):
                    out_rawFile.seek(headersize,1)
                    d1 = template_guppi_format(data[:,d1s*mb_i:d1s*(mb_i+1),:])
                    out_rawFile.write(d1.tostring())

            bend = time.time()
            print(f'block duration: {(bend-bstart)/60}')

        #===============================================
        #***********************************************

        #save output numpy arrays

        print(f'Flags: {self._flags_filename}')
        np.save(self._flags_filename, self.flags_all)
        print(f'Spect: {self._spect_filename}')
        np.save(self._spect_filename, self.spect_all)
        print(f'Regen: {self._regen_filename}')
        np.save(self._regen_filename, self.regen_all)


        if self.det_method == 'SK':
            print(f'SS-SK: {self._ss_sk_filename}')
            np.save(self._ss_sk_filename, self.ss_sk_all)
            print(f'MS-SK: {self._ms_sk_filename}')

            np.save(self._ms_sk_filename, self.ms_sk_all)
            #need to add the logging thing
            log = '/data/scratch/SKresults/SK_log.txt'
            os.system(f"""echo "'{self._spect_filename}','{self._flags_filename}','{self._regen_filename}','{self._ss_sk_filename}','{self._ms_sk_filename}'\n===============================" >> {log}""")

        elif self.det_method == 'AOF':
            #need to add the logging thing
            log = '/data/scratch/AOFresults/AOF_log.txt'
            os.system(f"""echo "'{self._spect_filename}','{self._flags_filename}','{self._regen_filename}'\n===============================" >> {log}""")

        elif self.det_method == 'MAD':
            log = '/data/scratch/MADresults/MAD_log.txt'
            os.system(f"""echo "'{self._spect_filename}','{self._flags_filename}','{self._regen_filename}'\n===============================" >> {log}""")

        elif self.det_method == 'SE':
            print(f'mod Z-score: {self._zsc_filename}')
            np.save(self._zsc_filename, self.zsc_all)
            log = '/data/scratch/SEresults/SE_log.txt'
            os.system(f"""echo "'{self._spect_filename}','{self._flags_filename}','{self._regen_filename}','{self._zsc_filename}'\n===============================" >> {log}""")
        

        # elif self.det_method == 'IQRM':
            
        #     print(f'avg pre: {self._avg_pre_filename}')
        #     np.save(self._avg_pre_filename, avg_pre)
            # print(f'avg post: {self._avg_post_filename}')
            # np.save(self._avg_post_filename, self.avg_post)


        #***********************************************
        #===============================================

        
        #flagging stuff
        template_print_flagstats(self.flags_all, True)



        #umm... done?
        end_time = time.time()
        dur = np.around((end_time-start_time)/60, 2)

        print(f'Duration: {dur} minutes')


    #resolution: output frequency resolution in kHz
    #mit: run on the mitigated or unmitigated data
    #mask: use the mask to skip over flagged bits when fine channelizing (not done)
#     def fine_channelize(self, resolution, mit=False, mask=False):
        
#         start_time = time.time()
#         if mit:
#             if mask:
#                 raw2spec_mask(resolution,self._outfile,mask)
#             else:
#                 raw2spec(resolution,self._outfile)
#         else:        
#             if mask:
#                 raw2spec_mask(resolution,self._rawFile,mask)
#             else:
#                 raw2spec(resolution,self._rawFile)
#         end_time = time.time()
#         dur = np.around((end_time-start_time)/60, 2)

#         print(f'Duration: {dur} minutes')

    def fine_channelize(self, resolution, mit=False, mask=False):
        start_time = time.time()
        if mask:
            if self.det_method=='AOF':
                in_mask = self.output_mit_srdp_dir+'*flags_block*.npy'
            else:
                in_mask = self._flags_filename
        if mit:
            if mask:
                out_fc_fname = f"{self.output_mit_srdp_dir}{self.npybase}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}_{resolution}_mask.spec.pkl"
                raw2spec_god(resolution,GuppiRaw(self.outfile_raw_full),self.det_method, out_fc_fname, in_mask)
            else:
                out_fc_fname = f"{self.output_mit_srdp_dir}{self.npybase}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}_{resolution}_nomask.spec.pkl"
                raw2spec_god(resolution,GuppiRaw(self.outfile_raw_full),self.det_method, out_fc_fname)
        else:
            if mask:
                out_fc_fname = f"{self.output_unmit_srdp_dir}{self.npybase}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}_{resolution}_mask.spec.pkl"
                raw2spec_god(resolution,self._rawFile,self.det_method, out_fc_fname, in_mask)
            else:
                out_fc_fname = f"{self.output_unmit_srdp_dir}{self.npybase}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}_{resolution}_nomask.spec.pkl"
                raw2spec_god(resolution,self._rawFile,self.det_method, out_fc_fname)
        end_time = time.time()
        dur = np.around((end_time-start_time)/60, 2)

        print(f'Duration: {dur} minutes')



    def pulsar_reduction(self):
        start_time = time.time()


        #reduce_pulsar_data(self.outfile_raw_full, self.output_srdp_dir)




        end_time = time.time()
        dur = np.around((end_time-start_time)/60, 2)

        print(f'Duration: {dur} minutes')



    def load_basic_srdps(self):
        

        # need to output: logs, (logsf or ave_f), logsm
        # for sk: sk, mssk

        #these are np.loads because you can do this function without first doing self.run_all()
        

        s = np.load(self._spect_filename)
        sm = np.load(self._regen_filename)

        logs = np.log10(s)
        logsm = np.log10(sm)

        out = (logs,logsm)
        
        if (self.det_method == 'AOF') or (self.det_method == 'MAD'):
            f = load_raw_flags(self.output_srdp_dir+self.npybase+'_flags_block*.npy',self.ave_factor)
            out = out + (f,)
        else:
            f = np.load(self._flags_filename)
            out = out + (f,)

        return out
            





























