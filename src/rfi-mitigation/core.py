
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

from utils import *

import iqrm

import RFI_detection as rfi
from tqdm import tqdm






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
    def __init__(self, det_method, infile, repl_method, cust, output_bool = True, mb=1, rawdata=False):
        #user-given attributes
        self.det_method = det_method        
        self.infile = template_infile_mod(infile,in_dir)[0]
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata  

        #default/hardcoded attributes
        self.in_dir = '/data/rfimit/unmitigated/rawdata/'#move to actual data dir
        self._rawFile = GuppiRaw(infile)


    def run_all(self):
        #do all the rfi mitigation steps
        #pass


    #def run_meta():
        #first step of RFI mitigation, get all the meta info set up correctly

       if self.output_bool:
	        template_check_outfile(infile,outfile)
	        out_rawFile = open(outfile,'rb+')

        
        template_check_nblocks(rawFile,mb)
        numblocks = self.rawFile.find_n_data_blocks()

        for block in range(numblocks//mb):
            print('------------------------------------------')
            print(f'Block: {(block*mb)+1}/{numblocks}')


            #print header for the first block
            if block == 0:
                template_print_header(rawFile)


            #loading multiple blocks at once?	
            for mb_i in range(mb):
                if mb_i==0:
                header,data = rawFile.read_next_data_block()
                data = np.copy(data)
                d1s = data.shape[1]
            else:
                h2,d2 = rawFile.read_next_data_block()
                data = np.append(data,np.copy(d2),axis=1)

            #find data shape
            num_coarsechan = data.shape[0]
            num_timesamples= data.shape[1]
            num_pol = data.shape[2]
            print(f'Data shape: {data.shape} || block size: {data.nbytes}')

            #save raw data?
            if rawdata:
                template_save_npy(data,block,npy_base)    
       

            #===============================================
            #***********************************************

            if self.method = 'sk':
                flags_block = sk-rfi.SK_detection(self,data)


        pass

    def run_mitigation():
        #heavy lifter, step through all blocks and do detection
        #replace and write back data to disk if given the true flag
        pass






























