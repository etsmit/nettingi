

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
        #user-given attributes
        self.det_method = 'IQRM'       
        self.infile = template_infile_mod(infile,self.in_dir)[0]
        self.repl_method = repl_method
        self.cust = cust
        # self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        # self.ave_factor = ave_factor 

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
            outfile_pattern = f'r{IQRM_radius}_t{IQRM_threshold}_{IQRM_datatype}'


        # any separate results filenames you need, in addition to the flags filename, put them here
        npybase = self.out_dir+'npy_results/'+infile[len(self.in_dir):-4]


        self._flags_filename = f"{npybase}_flags_{self.det_method}_{self._outfile_pattern}_{cust}.npy"

        self.avg_pre_filename = f"{npybase}_avg_pre_{self.det_method}_{outfile_pattern}_{cust}.npy"
        self.avg_post_filename = f"{npybase}_avg_post_{self.det_method}_{outfile_pattern}_{cust}.npy"
        self.spost_filename = f"{npybase}_spost_{self.det_method}_{outfile_pattern}_{cust}.npy"

        
    

        self._outfile = f"{self._jetstor_dir}{infile[len(self.in_dir):-4]}_{self.det_method}_{self._outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"
        #threshold calc from sigma
        self.IQRM_lag = iqrm.core.genlags(IQRM_radius, geofactor=1.5)
        print('integer lags, k: {}'.format(self.IQRM_lag))
        
        
