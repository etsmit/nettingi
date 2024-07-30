

import numpy as np
import os,sys
import matplotlib.pyplot as plt

import scipy as sp
import scipy.optimize
import scipy.special
import math as math

import argparse

import time

from blimpy import GuppiRaw


import iqrm

from tqdm import tqdm
from .core import mitigateRFI

from .utils import *
import iqrm

class rfi_iqrm(mitigateRFI):
    def __init__(self, infile, repl_method, IQRM_radius=5, IQRM_threshold=3.0, IQRM_datatype='std', IQRM_breakdown=512, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        valid = ["std", "power", "avg", "mad", "sk"] # valid inputs for IQRM_datatype
        if IQRM_datatype not in valid:
            raise ValueError("IQRM_datatype must be one of %r." % valid)
        
        #user-given attributes
        self.det_method = 'IQRM'  
        self.in_dir = '/data/rfimit/unmitigated/rawdata/'
        self.infile = template_infile_mod(infile,self.in_dir)[0]
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor 

        #default/hardcoded attributes
        self.in_dir = '/data/rfimit/unmitigated/rawdata/'#move to actual data dir
        self.infile = template_infile_mod(infile,self.in_dir)
        self._rawFile = GuppiRaw(self.infile)


        self.IQRM_radius = IQRM_radius
        self.IQRM_threshold = IQRM_threshold
        self.IQRM_datatype = IQRM_datatype
        self.IQRM_breakdown = IQRM_breakdown

        self._out_dir = '/data/scratch/IQRMresults/'
        self._jetstor_dir = '/jetstor/scratch/IQRM_rawdata_results/'


        if IQRM_datatype == 'std':
            self._outfile_pattern = f'r{IQRM_radius}_t{IQRM_threshold}_{IQRM_datatype}_b{IQRM_breakdown}'
        else:
            self._outfile_pattern = f'r{IQRM_radius}_t{IQRM_threshold}_{IQRM_datatype}'


        # any separate results filenames you need, in addition to the flags filename, put them here
        npybase = self._out_dir+'npy_results/'+infile[:-4]
        print(infile[:-4])
        


        self._flags_filename = f"{npybase}_flags_{self.det_method}_{self._outfile_pattern}_{cust}.npy"

        # self._avg_pre_filename = f"{npybase}_avg_pre_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        # self._avg_post_filename = f"{npybase}_avg_post_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        self._regen_filename = f"{npybase}_regen_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        self._spect_filename = f"{npybase}_spect_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        
    

        self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"
        #threshold calc from sigma
        IQRM_lag = iqrm.core.genlags(IQRM_radius, geofactor=1.5)
  
        print('integer lags, k: {}'.format(IQRM_lag))


        
    def iqrm_detection(self, data):
        """
        Performs the iqrm algorithm on the data. Transforms the data as necessary depending on the IQRM_datatype.
    
        Parameters
        -----------
        data : ndarray
            file name of the data that needs rfi mitigation
        Returns
        -----------
        out : ndarray
                An array of flags, indicating where the RFI is.
        """
        if self.IQRM_datatype == "std" and data.shape[1] % self.IQRM_breakdown != 0:
            raise ValueError("IQRM_breakdown must be a factor of %r." % self.shape[1])
            
        if data.shape[1] % self.ave_factor != 0:
            raise ValueError("ave_factor must be a factor of %r." % self.shape[1])
            
        if self.IQRM_datatype == 'power':
            flag_chunk = iqrm_power(data, self.IQRM_radius, self.IQRM_threshold)
    
        elif self.IQRM_datatype == 'avg':
            flag_chunk = iqrm_avg(data, self.IQRM_radius, self.IQRM_threshold, self.IQRM_breakdown)
            # standard dev
        else:# if self.IQRM_datatype == 'std': 
            flag_chunk = iqrm_std(data, self.IQRM_radius, self.IQRM_threshold, self.IQRM_breakdown)
    
        return flag_chunk
    
    
        # else:
        #     return #throw some error?
    
    # def averager(data,m):
    # 	"""
    # 	average of data
    # 	"""
    # 	step1_p0 = np.reshape(data[:,:,0], (data.shape[0],-1,m)) # polarization 1
    # 	step1_p1 = np.reshape(data[:,:,1], (data.shape[0],-1,m)) # 2
    # 	step2_p0 = np.expand_dims(np.mean(step1_p0,axis=2),axis=2)
    # 	step2_p1 = np.expand_dims(np.mean(step1_p1,axis=2),axis=2)
    # 	return np.concatenate((step2_p0,step2_p1),axis=2) 
                      
    
    # def stdever(data,m):
    # 	"""
    # 	standard deviation of data
    # 	"""
    # 	step1_p0 = np.reshape(data[:,:,0], (data.shape[0],-1,m)) # polarization 1
    # 	step1_p1 = np.reshape(data[:,:,1], (data.shape[0],-1,m)) # 2
    # 	step2_p0 = np.expand_dims(np.std(step1_p0,axis=2),axis=2)
    # 	step2_p1 = np.expand_dims(np.std(step1_p1,axis=2),axis=2)
    # 	return np.concatenate((step2_p0,step2_p1),axis=2)    
    
    
