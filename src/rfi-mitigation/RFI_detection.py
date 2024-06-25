
import numpy as np
import os,sys

import time

import matplotlib.pyplot as plt

from RFI_support import *

from numba import jit,prange
import iqrm
from tqdm import tqdm





#---------------------------------------------------------
# 1 . Functions for performing SK
#---------------------------------------------------------

#Compute SK on a 2D array of power values

@jit(parallel=True)
def SK_EST(a,m,n=1,d=1):
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

	#make s1 and s2 as defined by whiteboard (by 2010b Nita paper)
	a = a[:,:m]*n
	sum1=np.sum(a,axis=1)
	sum2=np.sum(a**2,axis=1)
	sk_est = ((m*n*d+1)/(m-1))*(((m*sum2)/(sum1**2))-1)                     
	return sk_est


@jit(parallel=True)
def SK_EST_alt(s1,s2,m,n=1,d=1):
	"""
	Compute SK on a 2D array of power values, using s1 and s2 given instead of data

	Parameters
	-----------

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
	sk_est = ((m*n*d+1)/(m-1))*(((m*s2)/(s1**2))-1)                     
	return sk_est


#multiscale variant
#only takes n=1 for now
#takes sum1 and sum2 as arguments rather than computing inside
@jit(parallel=True)
def ms_SK_EST(s1,s2,m,n=1,d=1):
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
	
                 #((m*n*d+1)/(m-1))*(((m*sum2)/(sum1**2))-1)
	sk_est = ((m*n*d+1)/(m-1))*(((m*s2)/(s1**2))-1)
	#print(sk_est)
	return sk_est


def master_SK_EST(data,SK_M,m,n=1,d=1):
	

	num_skbins = data.shape[1]//SK_M

	data = np.abs(np.reshape(data,(data.shape[0],-1,SK_M,data.shape[2])))**2


	s1 = np.sum(data,axis=2)
	s2 = np.sum(data**2,axis=2)

	spect_block = np.mean(data,axis=1)




def upperRoot(x, moment_2, moment_3, p):
	upper = np.abs( (1 - sp.special.gammainc( (4 * moment_2**3)/moment_3**2, (-(moment_3-2*moment_2**2)/moment_3 + x)/(moment_3/2/moment_2)))-p)
	return upper

#helps calculate lower SK threshold
def lowerRoot(x, moment_2, moment_3, p):
	lower = np.abs(sp.special.gammainc( (4 * moment_2**3)/moment_3**2, (-(moment_3-2*moment_2**2)/moment_3 + x)/(moment_3/2/moment_2))-p)
	return lower

#fully calculates upper and lower thresholds
#M = SK_ints
#default p = PFA = 0.0013499 corresponds to 3sigma excision
def SK_thresholds(M, N = 1, d = 1, p = 0.0013499):
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

	Nd = N * d
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
	upperThreshold = sp.optimize.newton(upperRoot, x[0], args = (moment_2, moment_3, p))
	lowerThreshold = sp.optimize.newton(lowerRoot, x[0], args = (moment_2, moment_3, p))
	return lowerThreshold, upperThreshold


	


def averager(data,m):
	"""
	average of data
	"""
	step1_p0 = np.reshape(data[:,:,0], (data.shape[0],-1,m)) # polarization 1
	step1_p1 = np.reshape(data[:,:,1], (data.shape[0],-1,m)) # 2
	step2_p0 = np.expand_dims(np.mean(step1_p0,axis=2),axis=2)
	step2_p1 = np.expand_dims(np.mean(step1_p1,axis=2),axis=2)
	return np.concatenate((step2_p0,step2_p1),axis=2)           
                  

def stdever(data,m):
	"""
	standard deviation of data
	"""
	step1_p0 = np.reshape(data[:,:,0], (data.shape[0],-1,m)) # polarization 1
	step1_p1 = np.reshape(data[:,:,1], (data.shape[0],-1,m)) # 2
	step2_p0 = np.expand_dims(np.mean(step1_p0,axis=2),axis=2)
	step2_p1 = np.expand_dims(np.mean(step1_p1,axis=2),axis=2)
	return np.concatenate((step2_p0,step2_p1),axis=2)    


def iqrm_power(data, radius, threshold):
	m = 512 # constant
	avg_pre = averager(np.abs(data)**2,m)
	data = np.abs(data)**2
	flag_chunk = np.zeros(data.shape)
	for i in tqdm(range(data.shape[2])): # iterate through polarizations
		for j in range(data.shape[0]): # iterate through channels
			flag_chunk[j,:,i] = iqrm.iqrm_mask(data[j,:,i], radius = radius, threshold = threshold)[0]
    
#     avg_post = 
	return flag_chunk, avg_pre



def iqrm_std(data, radius, threshold, breakdown):
	"""
	breakdown must be a factor of the time shape data[1].shape()
	"""
	m = 512 # constant
	avg_pre = averager(np.abs(data)**2, m)
# 	data_pol0 = stdever(np.abs(data[:,:,0])**2, breakdown) # make it a stdev
# # 	shape=np.expand_dims(shape, axis=2)
# 	flag_chunk = np.zeros((*data_pol0.shape[:2], 2))
# 	print('Data shape: {} || block size: {}'.format(flag_chunk.shape,flag_chunk.nbytes))
    
	data = stdever(np.abs(data)**2, breakdown)
	flag_chunk = np.zeros(data.shape)
	print('Flag shape: {} || block size: {}'.format(flag_chunk.shape,flag_chunk.nbytes))
	for i in tqdm(range(data.shape[2])): # iterate through polarizations
		for j in range(data.shape[0]): # iterate through channels
			flag_chunk[j,:,i] = iqrm.iqrm_mask(data[j,:,i], radius = radius, threshold = threshold)[0]
            
            
            
            
# 	for j in range(data_pol0.shape[0]): # iterate through channels
# 		flag_chunk[j,:,0] = iqrm.iqrm_mask(data_pol0[j,:], radius = radius, threshold = threshold)[0]
# 	data_pol1 = stdever(np.abs(data[:,:,1])**2, breakdown) # make it a stdev
# 	for j in range(data_pol1.shape[0]): # iterate through channels
# 		flag_chunk[j,:,1] = iqrm.iqrm_mask(data_pol1[j,:], radius = radius, threshold = threshold)[0]

	return flag_chunk, avg_pre



















