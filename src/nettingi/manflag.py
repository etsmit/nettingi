#h
import numpy as np

import math as math
from blimpy import GuppiRaw
import matplotlib.pyplot as plt
import time
from .plotting import implot



from .utils import (
    template_bookkeeping,
    template_check_outfile,
    template_calc_ave,
    template_check_nblocks,
    template_print_header,
    )


class rfi_manual():
    #h
    def __init__(self, infile, repl_method, cust='', output_bool = True, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'MAN'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = 1
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile 



        self._outfile_pattern = f"{self.repl_method}_{self.cust}"    
        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir, self.output_unmit_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern,self.det_method)
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]
        

        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"


        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)


    def manualflag(self):

        start_time = time.time()
        if self.output_bool:
            template_check_outfile(self.infile_raw_full,self.outfile_raw_full)
            out_rawFile = open(self.outfile_raw_full,'rb+')

        
        template_check_nblocks(self._rawFile,self.mb)
        numblocks = self._rawFile.find_n_data_blocks()

        for bi in range(numblocks//self.mb):
            print('------------------------------------------')
            print(f'Block: {(bi*self.mb)+1}/{numblocks}')


            #print header for the first block
            if bi == 0:
                headersize = template_print_header(self._rawFile)


            #loading multiple blocks at once?        
            for mb_i in range(self.mb):
                if mb_i==0:
                    header,data = self._rawFile.read_next_data_block()
                    self.bw = int(header["OBSBW"])
                    self.rf = float(header["OBSFREQ"])
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

            if bi == 0:
                self.spect_all = spect_block
            else:
                self.spect_all = np.concatenate((self.spect_all, spect_block),axis=1)

        
        self.logs = np.log10(self.spect_all)

        plot_and_flag()


    def plot_and_flag(self):

            ax_lw = 3
            ms = s.shape[1]*(self.ave_factor*1e3*(s.shape[0]/(bw*1e6)))
            ext = [0,ms/1e3,rf-0.5*bw,rf+0.5*bw]

            self._fig, self._ax = plt.subplots(figsize=(14,10))
            ax_lw = 3
            self._ax.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2)
            self._ax.spines['bottom'].set_linewidth(ax_lw)
            self._ax.spines['top'].set_linewidth(ax_lw)
            self._ax.spines['left'].set_linewidth(ax_lw)
            self._ax.spines['right'].set_linewidth(ax_lw)

            self._ax.imshow(s,interpolation='nearest',aspect='auto',cmap='hot',vmin=vmin,vmax=vmax,extent=ext)
            plt.ylabel('Frequency (MHz)',fontsize=20)
            plt.xlabel('Time (sec)',fontsize=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            cbar = plt.colorbar()
            for t in cbar.ax.get_yticklabels():
                t.set_fontsize(20)
            plt.tight_layout()
            plt.show()

        #implot(self.logs[:,:,0])









