
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
    def __init__(self, det_method, infile, repl_method, cust):
        #user-given attributes
        self.det_method = det_method        
        self.infile = infile
        self.repl_method = repl_method
        self.cust = cust

    
        #default/hardcoded attributes (these should be changed)
        self.in_dir = '/data/rfimit/unmitigated/rawdata/'
        self.out_dir = '/data/scratch/SKresults/'
        self.jetstor_dir = '/jetstor/scratch/SK_rawdata_results/'




    def run_all():
        #do all the rfi mitigation steps
        pass


    def run_meta():
        #first step of RFI mitigation, get all the meta info set up correctly
        pass

    def run_mitigation():
        #heavy lifter, step through all blocks and do detection
        #replace and write back data to disk if given the true flag
        pass






























