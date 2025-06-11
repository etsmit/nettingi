#h
import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
import torch

import nettingi

from blimpy import GuppiRaw

from .core import mitigateRFI
from numba import prange
from ConvRFI  import  RFIconv,init_RFIconv

from .utils import (
    template_bookkeeping,
    )
from .utils import template_calc_ave
#for object-based JIT compile
# spec = [
#  'SK_detection',
#  'SK_m',
#  'SK_thresholds',
#  '_SK_p',
#  '_flags_filename',
#  '_jetstor_dir',
#  '_lt',
#  '_ms_lt',
#  '_ms_sk_filename',
#  '_ms_ut',
#  '_out_dir',
#  '_outfile',
#  '_outfile_pattern',
#  '_rawFile',
#  '_regen_filename',
#  '_spect_filename',
#  '_ss_sk_filename',
#  '_ut',
#  'ave_factor',
#  'cust',
#  'd',
#  'det_method',
#  'flags_all',
#  'in_dir',
#  'infile',
#  'lowerRoot',
#  'mb',
#  'ms0',
#  'ms1',
#  'ms_sk_all',
#  'mssk',
#  'multi_scale_SK_EST',
#  'n',
#  'output_bool',
#  'rawdata',
#  'regen_all',
#  'repl_method',
#  'run_all',
#  'sigma',
#  'single_scale_SK_EST',
#  'spect_all',
#  'ss_sk_all',
#  'upperRoot']

class rfi_ConvRFI(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, a0, a1, a2, a3, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'Conv'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile 


        #convRFI related parameters
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        

        self._outfile_pattern = f"a0{self.a0}_a1{self.a1}_a2{self.a2}-_a3{self.a3}_{self.repl_method}_{self.cust}"    
        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir, self.output_unmit_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern,self.det_method)
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]
        

        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        self._ss_sk_filename = f"{self.output_mit_srdp_dir}{self.npybase}_skval_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._ms_sk_filename = f"{self.output_mit_srdp_dir}{self.npybase}_mssk_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"


        #self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"
        

        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)


#@jit(nopython=True, parallel=True)
    def conv_detection(self, data):
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        if device!='cpu':
            torch.cuda.empty_cache()
        
        net = RFIconv(device=device)
        agg_factor = [self.a0,self.a1,self.a2,self.a3]

        s_ave = template_calc_ave(data,self.ave_factor).astype(np.float32)
        net = init_RFIconv(net,aggressive_factor = agg_factor,device=device).to(device)
        flags_block = np.empty(s_ave.shape)
        for pol in range(s_ave.shape[2]):
            with torch.no_grad():
                output = net(torch.tensor(s_ave[:,:,pol].squeeze()[None,None,:,:]).to(device)).squeeze().cpu().numpy()
            flags_block[:,:,pol] = output
            
        net = init_RFIconv(net,aggressive_factor = agg_factor,device=device).to(device)
        
        if device!='cpu':
            torch.cuda.empty_cache() 
        return flags_block 