#spectral entropy

import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
from blimpy import GuppiRaw
import matplotlib.pyplot as plt

from .core import mitigateRFI

from .utils import *

from numba import jit

from tqdm import tqdm


class rfi_se(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, m, s, nbits=8, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'SE'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile 

        #default/hardcoded attributes

        self.sigma = s
        self.SE_m = m
        self.nbits = nbits

        self._outfile_pattern = f"m{self.SE_m}_s{self.sigma}"    

        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir,self.output_unmit_srdp_dir = template_bookkeeping(self.infile,self._outfile_pattern,self.det_method)
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]
        
        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"


        #self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"
        
        
        #any derived thresholds/arrays
        self._zsc_filename = f"{self.output_mit_srdp_dir}{self.npybase}_zsc_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)



    def se_detection(self,data):
        
        nbit_range = np.arange(2**self.nbits/-2 + 1, 2**self.nbits/2) 

        #get H for both real and complex

        a = np.reshape(data,(data.shape[0],-1,self.SE_m,data.shape[2]))

        H = np.empty((a.shape[0],a.shape[1],a.shape[3]),dtype=np.complex64)
        Hbase = np.empty((a.shape[0],a.shape[1],a.shape[3]),dtype=np.complex64)


        print('starting H...')
        for pol in range(a.shape[3]):
            for chan in tqdm(range(a.shape[0])):
                for tb in range(a.shape[1]):
                    H[chan,tb,pol], Hbase[chan,tb,pol] = getHs(a[chan,tb,:,pol], nbit_range)
        print('sre')
        #calculate relative spectral entropy
        sre = np.abs(H - Hbase)

        print('zscore...')
        #calculate significance of outliers
        modZscore = get_modZscore(sre)
        
        print('flagging')
        #flag
        flags_block = np.zeros(sre.shape,dtype=np.int8)
        flags_block[np.isnan(sre)] = 1
        flags_block[modZscore > self.sigma] = 1

        return flags_block, modZscore



def getHs(a,bins):
    numr,bins,patches = plt.hist(a.real, bins=bins)
    numr /= np.sum(numr)
    pdfr = get_gaussfit(bins,numr)
    numi,bins,patches = plt.hist(a.imag, bins=bins)
    numi /= np.sum(numi)
    pdfi = get_gaussfit(bins,numi)
    outH = np.sum(-pdfr*np.log(pdfr)) + 1.j*np.sum(-pdfi*np.log(pdfi))
    outHbase=(0.5*(1+np.log(2*np.pi*np.var(a.real)))) + 1.j*(0.5*(1+np.log(2*np.pi*np.var(a.imag))))

    return outH, outHbase

def gauss_func(x,a,b):
    return a*np.exp(-b*x**2)


def get_gaussfit(bins,num):
    #need to return nans in case fit fails. will flag these pixels back in parent function
    try:
        popt,pcov = scipy.optimize.curve_fit(gauss_func,bins[:-1],num,p0=(0,1))
        return gauss_func(bins[:-1],popt[0],popt[1])
    except:
        popt = np.array([np.nan,np.nan])
        return popt

def get_modZscore(a):

    mad = np.median(np.abs(a-np.median(a,axis=1)))
    M = (0.6745*(a-np.median(a,axis=1)))/mad

    return M










