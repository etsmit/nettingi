

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
    def __init__(self, infile, repl_method, IQRM_radius=5, IQRM_threshold=3.0, IQRM_datatype='std', IQRM_breakdown=512, rawdata=False):
        #user-given attributes
        self.det_method = 'IQRM'       
        self.infile = template_infile_mod(infile,self.in_dir)[0]
        self.repl_method = repl_method
        # self.cust = cust
        # self.output_bool = output_bool 
        # self.mb = mb
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


        self._flags_filename = f"{npybase}_flags_{IDstr}_{outfile_pattern}_{cust}.npy"

        self._sk_filename = f"{npybase}_skval_{IDstr}_{outfile_pattern}_{cust}.npy"
        self._mssk_filename = f"{npybase}_mssk_{IDstr}_{outfile_pattern}_{cust}.npy"


        self._outfile = f"{jetstor_dir}{infile[len(in_dir):-4]}_{IDstr}_{outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"
        
        

        #any derived thresholds/arrays
        self._SK_p = (1-scipy.special.erf(sigma/math.sqrt(2))) / 2
        print(f'Probability of false alarm: {SK_p}')

        #calculate thresholds
        print('Calculating SK thresholds...')
        self._lt, self._ut = self.SK_thresholds(self)
        print(f'Upper Threshold: {ut}'+str(ut))
        print(f'Lower Threshold: {lt}'+str(lt))

        #calculate ms thresholds
        self._ms_lt, self._ms_ut = self.SK_thresholds(SK_M*ms0*ms1, N = n, d = d, p = SK_p)
        print(f'MS Upper Threshold: {ms_ut}')
        print(f'MS Lower Threshold: {ms_lt}')