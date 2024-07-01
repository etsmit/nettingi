#h



import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
from blimpy import GuppiRaw

from .core import mitigateRFI

from .utils import *


class rfi_sk(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, m, mssk, n, d, s, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'SK'       
        #self.infile = template_infile_mod(infile,self.in_dir)[0]
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


        self.SK_m = m
        self.mssk = mssk
        self.n = n
        self.d = d
        self.sigma = s

        self._out_dir = '/data/scratch/SKresults/'
        self._jetstor_dir = '/jetstor/scratch/SK_rawdata_results/'

        self.ms0 = int(mssk.split(',')[0])
        self.ms1 = int(mssk.split(',')[1])

        
        self._outfile_pattern = f"m{self.SK_m}_s{self.sigma}_ms{self.ms0}-{self.ms1}"


        # any separate results filenames you need, in addition to the flags filename, put them here
        npybase = self._out_dir+'npy_results/'+infile[len(self.in_dir):-4]


        self._flags_filename = f"{npybase}_flags_{self.det_method}_{self._outfile_pattern}_{self.cust}.npy"

        self._sk_filename = f"{npybase}_skval_{self.det_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._mssk_filename = f"{npybase}_mssk_{self.det_method}_{self._outfile_pattern}_{self.cust}.npy"


        self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"
        
        

        #any derived thresholds/arrays
        self._SK_p = (1-scipy.special.erf(self.sigma/math.sqrt(2))) / 2
        print(f'Probability of false alarm: {self._SK_p}')

        #calculate thresholds
        print('Calculating SK thresholds...')
        self._lt, self._ut = self.SK_thresholds(self.SK_m)
        print(f'Upper Threshold: {self._ut}')
        print(f'Lower Threshold: {self._lt}')

        #calculate ms thresholds
        self._ms_lt, self._ms_ut = self.SK_thresholds(self.SK_m*self.ms0*self.ms1)
        print(f'MS Upper Threshold: {self._ms_ut}')
        print(f'MS Lower Threshold: {self._ms_lt}')


    def SK_detection(self,data):

        s = np.abs(data)**2

        num_coarsechan = s.shape[0]
        num_timesamples = s.shape[1]
        num_SKbins = num_timesamples // self.SK_m

        flags_block = np.zeros((s.shape[0], s.shape[1]//self.SK_m, s.shape[2]),dtype=np.int8)
        ms_flags_block = np.zeros((s.shape[0] - (self.ms0-1) , flags_block.shape[1] - (self.ms1-1) , s.shape[2]))

        #detect and flag RFI
        for pol in range(s.shape[2]):


            #single scale
            ss_sk_block = self.single_scale_SK_EST(s)
            flags_block[ss_sk_block < self._lt] = 1
            flags_block[ss_sk_block > self._ut] = 1


            #multiscale
            if (self.mssk != '0,0'):
                ms_sk_block = self.multi_scale_SK_EST(s)
                ms_flags_block[ms_sk_block < self._ms_lt] = 1
                ms_flags_block[ms_sk_block > self._ms_ut] = 1


                for ichan in range(self.ms0):
                    for itime in range(self.ms1):
                        flags_block[ichan:ichan+(num_coarsechan-(self.ms0-1)),itime:itime+(num_SKbins-(self.ms1-1)),:][ms_flags_block==1] = 1

        #track intermediate results
        #if bi == 0:
        #    self.flags_all = flags_block
        #    self.ss_sk_all = ss_sk_block
        #    self.ms_sk_all = ms_sk_block
        #else:
        #    self.ss_sk_all = np.concatenate((self.ss_sk_all, flags_block),axis=1)
        #    self.ss_sk_all = np.concatenate((self.ss_sk_all, ss_sk_block),axis=1)
        #    self.ms_sk_all = np.concatenate((self.ms_sk_all, ms_sk_block),axis=1)



        return flags_block,ss_sk_block,ms_sk_block
        
        #pass



    def single_scale_SK_EST(self,s):
        """
        Compute SK on a 2D array of power values.

        Parameters
        -----------
        a : ndarray
            2-dimensional array of power values. Shape (Num Channels , Num Raw Spectra)
        n : int
            integer value of N in the SK function. Inside accumulations of spectra.
        m : int
            integer value of M in the SK function. Outside accumulations of spectra.
        d : float
            shape parameter d in the SK function. Usually 1 but can be empirically determined.
    
        Returns
        -----------
        out : ndarray
            Spectrum of SK values.
        """
        nd = self.n * self.d
        a = np.reshape(s,(s.shape[0],-1,self.SK_m,s.shape[2]))
        sum1=np.sum(a,axis=2)
        sum2=np.sum(a**2,axis=2)
        sk_est = ((self.SK_m*nd+1)/(self.SK_m-1))*(((self.SK_m*sum2)/(sum1**2))-1)                     
        return sk_est



    def multi_scale_SK_EST(self,s):
        """
        Multi-scale Variant of SK_EST.

        Parameters
        -----------
        s1 : ndarray
                2-dimensional array of power values. Shape (Num Channels , Num Raw Spectra)

        s2 : ndarray
                2-dimensional array of squared power values. Shape (Num Channels , Num Raw Spectra)

        m : int
                integer value of M in the SK function. Outside accumulations of spectra.

        ms0 : int
                axis 0 multiscale
        
        ms1 : int
                axis 1 multiscale
        
        Returns
        -----------
        out : ndarray
                Spectrum of SK values.
        """

        ms_binsize = self.ms0*self.ms1
        num_SKbins = s.shape[1] // self.SK_m
        ms_s1 = np.zeros((s.shape[0]-(self.ms0-1),num_SKbins-(self.ms1-1),2))
        ms_s2 = np.zeros((s.shape[0]-(self.ms0-1),num_SKbins-(self.ms1-1),2))
        nd = self.n * self.d

        a = np.reshape(s,(s.shape[0],-1,self.SK_m,s.shape[2]))
        s1=np.sum(a,axis=2)
        s2=np.sum(a**2,axis=2)


        #make multiscale S1, S2
        for ichan in range(self.ms0):
                for itime in range(self.ms1):
                        ms_s1 += (1./ms_binsize) * (s1[ichan:ichan+(s.shape[0]-(self.ms0-1)),itime:itime+(num_SKbins-(self.ms1-1)),:])
                        ms_s2 += (1./ms_binsize) * (s2[ichan:ichan+(s.shape[0]-(self.ms0-1)),itime:itime+(num_SKbins-(self.ms1-1)),:])
        
                 #((m*n*d+1)/(m-1))*(((m*sum2)/(sum1**2))-1)




        sk_est = ((self.SK_m*nd+1)/(self.SK_m-1))*(((self.SK_m*ms_s2)/(ms_s1**2))-1)
        #print(sk_est)
        return sk_est




    def upperRoot(self,x, moment_2, moment_3, p):
        upper = np.abs( (1 - sp.special.gammainc( (4 * moment_2**3)/moment_3**2, (-(moment_3-2*moment_2**2)/moment_3 + x)/(moment_3/2/moment_2)))-p)
        return upper

    #helps calculate lower SK threshold
    def lowerRoot(self,x, moment_2, moment_3, p):
        lower = np.abs(sp.special.gammainc( (4 * moment_2**3)/moment_3**2, (-(moment_3-2*moment_2**2)/moment_3 + x)/(moment_3/2/moment_2))-p)
        return lower

    #fully calculates upper and lower thresholds
    #M = SK_ints
    #default p = PFA = 0.0013499 corresponds to 3sigma excision
    def SK_thresholds(self,M):
        """
        Determine SK thresholds numerically.

        Parameters
        -----------
        m : int
                integer value of M in the SK function. Outside accumulations of spectra.
        n : int
                integer value of N in the SK function. Inside accumulations of spectra.
        d : float
                shape parameter d in the SK function. Usually 1 but can be empirically determined.
        p : float
                Prob of false alarm. 0.0013499 corresponds to 3-sigma excision.
        
        Returns
        -----------
        out : tuple
                Tuple of (lower threshold, upper threshold).
        """

        Nd = self.n * self.d
        p = self._SK_p
        #Statistical moments
        moment_1 = 1
        moment_2 = float(( 2*(M**2) * Nd * (1 + Nd) )) / ( (M - 1) * (6 + 5*M*Nd + (M**2)*(Nd**2)) )
        moment_3 = float(( 8*(M**3)*Nd * (1 + Nd) * (-2 + Nd * (-5 + M * (4+Nd))) )) / ( ((M-1)**2) * (2+M*Nd) *(3+M*Nd)*(4+M*Nd)*(5+M*Nd))
        moment_4 = float(( 12*(M**4)*Nd*(1+Nd)*(24+Nd*(48+84*Nd+M*(-32+Nd*(-245-93*Nd+M*(125+Nd*(68+M+(3+M)*Nd)))))) )) / ( ((M-1)**3)*(2+M*Nd)*(3+M*Nd)*(4+M*Nd)*(5+M*Nd)*(6+M*Nd)*(7+M*Nd) )
        #Pearson Type III Parameters
        delta = moment_1 - ( (2*(moment_2**2))/moment_3 )
        beta = 4 * ( (moment_2**3)/(moment_3**2) )
        alpha = moment_3 / (2 * moment_2)
        beta_one = (moment_3**2)/(moment_2**3)
        beta_two = (moment_4)/(moment_2**2)
        error_4 = np.abs( (100 * 3 * beta * (2+beta) * (alpha**4)) / (moment_4 - 1) )
        kappa = float( beta_one*(beta_two+3)**2 ) / ( 4*(4*beta_two-3*beta_one)*(2*beta_two-3*beta_one-6) )
        print('kappa: {}'.format(kappa))
        x = [1]
        print(x, moment_2, moment_3, p)
        upperThreshold = sp.optimize.newton(self.upperRoot, x[0], args = (moment_2, moment_3, p))
        lowerThreshold = sp.optimize.newton(self.lowerRoot, x[0], args = (moment_2, moment_3, p))
        return lowerThreshold, upperThreshold





        
