#mad

import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
from blimpy import GuppiRaw

from .core import mitigateRFI

from .utils import *

from numba import jit




class rfi_mad(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, m, s, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
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

        self._outfile_pattern = f"m{self.mad_m}_s{self.sigma}"    

        self.infile_raw_full, self.outfile_raw_full, self.output_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern)
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]
        
        self._flags_filename = f"{self.output_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"


        #self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"
        
        
        #any derived thresholds/arrays


        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)




    def mad_detection(self,data):



        s = np.abs(data)**2
        a = np.reshape(s,(s.shape[0],-1,self.MAD_m,s.shape[2]))

        median = np.median(a,axis=2)
        mad = np.median(np.abs(a-median))

        pulse = np.ones((1,self.MAD_m,1))

        sigma_r = np.kron( 1.4826 * mad, pulse )
        

        f = np.zeros(s.shape,dtype=np.int8)
        

    





