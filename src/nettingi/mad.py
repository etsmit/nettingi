#mad

import numpy as np

import math as math
from blimpy import GuppiRaw

from .core import mitigateRFI

from .utils import template_bookkeeping


class rfi_mad(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, m=512, n=128, s=3.0, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'MAD'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile 

        #default/hardcoded attributes

        self.sigma = s
        self.MAD_m = m
        self.MAD_n = n

        self._outfile_pattern = f"m{self.MAD_m}_n{self.MAD_n}_s{3.0}"    

        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern,self.det_method)
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]
        
        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        self._ut_filename = f"{self.output_mit_srdp_dir}{self.npybase}_ut_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._lt_filename = f"{self.output_mit_srdp_dir}{self.npybase}_lt_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"



        #self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"
        
        
        #any derived thresholds/arrays


        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)




    def mad_detection(self,data):



        s = np.abs(data)**2
        a = np.reshape(s,(s.shape[0],-1,self.MAD_m,s.shape[2]))

        pulse = np.ones((1,self.MAD_m,1))

        median = np.kron( np.median(a,axis=2), pulse )
        mad = np.median(np.abs(s-median))

        sigma_r = 1.4826 * mad
        ut = median + sigma_r
        lt = median - sigma_r

        f = np.zeros(s.shape,dtype=np.int8)
        
        f[s > ut] = 1
        f[s < lt] = 1


        return f


    def mad_detection_inside(self,data,N,M,th_mod):

        out_shape = (data.shape[0],data.shape[1]//N,data.shape[2])
        out_f = np.zeros(out_shape)
        out_ut = np.zeros(out_shape)
        out_lt = np.zeros(out_shape)

        for i in range(data.shape[2]):

            td = data[:,:,i]

            if td.shape[1] // (N*M) != td.shape[1] / (N*M):
                print(f'{N} x {M} need to integer divide {td.shape[1]}')
                exit()

            
            s = np.abs(td)**2
            a = np.mean(np.reshape(s,(s.shape[0],-1,N)),axis=2)

            b = np.reshape(a,(a.shape[0],-1,M))

            Mpulse = np.ones((1,1,M))
            # Npulse = np.ones((1,N))

            median = np.kron( np.expand_dims(np.median(b,axis=2),axis=2), Mpulse )
            #median = np.kron( np.median(b,axis=2), Mpulse )
            mad = np.kron( np.expand_dims(np.median(np.abs(b-median),axis=2),axis=2), Mpulse )

            sigma_r = (1.4826*th_mod) * mad
            ut = median + sigma_r
            lt = median - sigma_r

            f = np.zeros(b.shape,dtype=np.int8)
                
            f[b > ut] = 1
            f[b < lt] = 1

            f = np.reshape(f,(f.shape[0],f.shape[1]*M))
            #f = np.kron( f, Npulse)

            ut = np.reshape(ut,(ut.shape[0],ut.shape[1]*M))
            #ut = np.kron( ut, Npulse)

            lt = np.reshape(lt,(lt.shape[0],lt.shape[1]*M))
            #lt = np.kron( lt, Npulse)

            out_f[:,:,i] = f
            out_ut[:,:,i] = ut
            out_lt[:,:,i] = lt

            #print(f'f shape: {f.shape}')

        return out_f,out_ut,out_lt




