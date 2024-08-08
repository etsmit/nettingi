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
    def __init__(self, infile, repl_method, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        # valid = ["std", "power", "avg", "mad", "sk"] # valid inputs for IQRM_datatype
        # if IQRM_datatype not in valid:
        #     raise ValueError("IQRM_datatype must be one of %r." % valid)
        
        #user-given attributes
        self.det_method = 'AOF'  
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


        # self.IQRM_radius = IQRM_radius
        # self.IQRM_threshold = IQRM_threshold
        # self.IQRM_datatype = IQRM_datatype
        # self.IQRM_breakdown = IQRM_breakdown

        self._out_dir = '/data/scratch/SKresults/'
        self._jetstor_dir = '/jetstor/scratch/SK_rawdata_results/'


        # if IQRM_datatype == 'std':
        #     self._outfile_pattern = f'r{IQRM_radius}_t{IQRM_threshold}_{IQRM_datatype}_b{IQRM_breakdown}'
        # else:
        self._outfile_pattern = f''


        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self._out_dir+'npy_results/'+infile[:-4]


        self._flags_filename = f"{self.npybase}_flags_{self.det_method}_{self._outfile_pattern}_{cust}.npy"

        # self._avg_pre_filename = f"{npybase}_avg_pre_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        # self._avg_post_filename = f"{npybase}_avg_post_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        self._regen_filename = f"{self.npybase}_regen_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        self._spect_filename = f"{self.npybase}_spect_{self.det_method}_{self._outfile_pattern}_{cust}.npy"
        
    

        self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"
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
        flag_chunk = aof(data)
    
        return flag_chunk



def aof(data):
    nch = data.shape[0]
    ntimes = data.shape[1]//1
    count = 2
    
    aoflag = aoflagger.AOFlagger()
    
    # Load strategy from disk (alternatively use 'make_strategy' to use a default one)
    path = aoflag.find_strategy_file(aoflagger.TelescopeId.Generic)
    strategy = aoflag.load_strategy_file(path)
    
    aof_data = aoflag.make_image_set(ntimes, nch, count)
    
    print("Number of times: " + str(aof_data.width()))
    print("Number of channels: " + str(aof_data.height()))
    
    # When flagging multiple baselines, iterate over the baselines and
    # call the following code for each baseline
    # (to use multithreading, make sure to create an imageset for each
    # thread)
    flags_block = np.zeros(data.shape,dtype = np.int8)
    
    # Divide up the block into 32 time chunks for lighter RAM usage
    tb_size = data.shape[1]//1    

    for tb in tqdm(range(1)):
        tstart = tb*tb_size
        tend = (tb+1)*tb_size
        # Make 4 images: real and imaginary for2 pol
        for pol in range(data.shape[2]):
            aof_data.set_image_buffer(0,(data[:,tstart:tend,pol].real).astype(np.int8))
            aof_data.set_image_buffer(1, (data[:,tstart:tend,pol].imag).astype(np.int8))
    
            flags = strategy.run(aof_data)
        
        # flagvalues = flags.get_buffer()
        # flagcount = sum(sum(flagvalues))
        # print(
        #     "Percentage flags on zero data: "
        #     + str(flagcount * 100.0 / (nch * ntimes))
        #     + "%"
        # )
            flags_block[:,tstart:tend,pol] = flags.get_buffer()
    # flags.x = flags
    # flags = id(flags)
    # print(flags)
    # Collect statistics
    # We create some unrealistic time and frequency arrays to be able
    # to run these functions. Normally, these should hold the time
    # and frequency values.


    # flagger = aoflagger.AOFlagger()
    # path = flagger.find_strategy_file(aoflagger.TelescopeId.Generic)
    # strategy = flagger.load_strategy_file(path)
    # data1 = flagger.make_image_set(ntimes, nch, 8)

    # aoflagger.FlagMask()

    
    # ratiosum = 0.0
    # ratiosumsq = 0.0
    # for repeat in range(count):
    #     for imgindex in range(8):
    #         # Initialize data with random numbers
    #         values = data
    #         data1.set_image_buffer(imgindex, values)
    
    #     flags = strategy.run(data)
    #     flagvalues = flags.get_buffer()
    #     ratio = float(sum(sum(flagvalues))) / (nch*ntimes)
    #     ratiosum += ratio
    #     ratiosumsq += ratio*ratio
    
    # print("Percentage flags (false-positive rate) on Gaussian data: " +
    #     str(ratiosum * 100.0 / count) + "% +/- " +
    #     str(np.sqrt(
    #         (ratiosumsq/count - ratiosum*ratiosum / (count*count) )
    #         ) * 100.0) )
    return flags_block








