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


import aoflagger

from tqdm import tqdm
from .core import mitigateRFI

from .utils import *

class rfi_aof(mitigateRFI):
    def __init__(self, infile, strategy, num_images, repl_method, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        # valid = ["std", "power", "avg", "mad", "sk"] # valid inputs for IQRM_datatype
        # if IQRM_datatype not in valid:
        #     raise ValueError("IQRM_datatype must be one of %r." % valid)
        
        #user-given attributes
        self.det_method = 'AOF'  
        self.in_dir = '/data/rfimit/unmitigated/rawdata/'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile
        self.strategy = strategy
        self.num_images = int(num_images)

        #default/hardcoded attributes
        #self.in_dir = '/data/rfimit/unmitigated/rawdata/'#move to actual data dir
        #self.infile = template_infile_mod(infile,self.in_dir)
        #self._rawFile = GuppiRaw(self.infile)


        # self.IQRM_radius = IQRM_radius
        # self.IQRM_threshold = IQRM_threshold
        # self.IQRM_datatype = IQRM_datatype
        # self.IQRM_breakdown = IQRM_breakdown

        self._out_dir = '/data/scratch/SKresults/'
        self._jetstor_dir = '/jetstor/scratch/SK_rawdata_results/'


        # if IQRM_datatype == 'std':
        #     self._outfile_pattern = f'r{IQRM_radius}_t{IQRM_threshold}_{IQRM_datatype}_b{IQRM_breakdown}'
        # else:
        self._outfile_pattern = f'{self.strategy}_{self.num_images}images'
        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir, self.output_unmit_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern,self.det_method)
        self._rawFile = GuppiRaw(self.infile_raw_full)


        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]


        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        
    

        #self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"
        #threshold calc from sigma


        
    def aof_detection(self, data):
        """
        Performs the aof algorithm on the data. Transforms the data as necessary depending on the IQRM_datatype.
    
        Parameters
        -----------
        data : ndarray
            file name of the data that needs rfi mitigation
        Returns
        -----------
        out : ndarray
                An array of flags, indicating where the RFI is.
        """

            
        # if self.IQRM_datatype == 'power':
        #     flag_chunk = iqrm_power(data, self.IQRM_radius, self.IQRM_threshold)
    
        # elif self.IQRM_datatype == 'avg':
        #     flag_chunk = iqrm_avg(data, self.IQRM_radius, self.IQRM_threshold, self.IQRM_breakdown)
        #     # standard dev
        # else:# if self.IQRM_datatype == 'std': 
        flag_chunk = aof(data,self.num_images,self.strategy)
    
        return flag_chunk



def aof(data,count,strategy):
    nch = data.shape[0]
    ntimes = data.shape[1]
    #count = self.num_images
    
    aoflag = aoflagger.AOFlagger()
    
    # Load strategy from disk (alternatively use 'make_strategy' to use a default one)
    if strategy == 'generic':
        path = aoflag.find_strategy_file(aoflagger.TelescopeId.Generic)
    elif strategy == 'arecibo':
        path = aoflag.find_strategy_file(aoflagger.TelescopeId.Arecibo)
    elif strategy == 'parkes':
        path = aoflag.find_strategy_file(aoflagger.TelescopeId.Parkes)
    #print(path)
    strategy = aoflag.load_strategy_file(path)
    
    
    #print("Number of times: " + str(aof_data.width()))
    #print("Number of channels: " + str(aof_data.height()))
    
    # When flagging multiple baselines, iterate over the baselines and
    # call the following code for each baseline
    # (to use multithreading, make sure to create an imageset for each
    # thread)
    flags_block = np.zeros(data.shape,dtype = np.int8)
    
    if count==2:

        aof_data = aoflag.make_image_set(ntimes, nch, count)


        for pol in range(data.shape[2]):
            aof_data.set_image_buffer(0,(data[:,:,pol].real).astype(np.int8))
            aof_data.set_image_buffer(1,(data[:,:,pol].imag).astype(np.int8))
            #aof_data = make_aof_data(data[:,:,pol],aoflag,ntimes,nch,count)
    
            flags_block[:,:,pol] = (strategy.run(aof_data)).get_buffer()

    else:
        aof_data = make_aof_data(data,aoflag,ntimes,nch,count)
    
        flags_block = strategy.run(aof_data)
        

    #flags_block = flags.get_buffer()
    # flags.x = flags

    return flags_block

def make_aof_data(data,aoflag,ntimes,nch,count):

    aof_data = aoflag.make_image_set(ntimes, nch, count)

    if count==1:
        aof_data.set_image_buffer(0,np.abs(data))
        
    elif count==2:
        aof_data.set_image_buffer(0,(data.real).astype(np.int8))
        aof_data.set_image_buffer(1, (data.imag).astype(np.int8))


    elif count==4:
        
        aof_data.set_image_buffer(0,(data[:,:,0].real).astype(np.int8))
        aof_data.set_image_buffer(1, (data[:,:,0].imag).astype(np.int8))
        aof_data.set_image_buffer(2,(data[:,:,1].real).astype(np.int8))
        aof_data.set_image_buffer(3, (data[:,:,1].imag).astype(np.int8))
















