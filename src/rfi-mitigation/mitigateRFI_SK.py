"""
mitigateRFI


Template for RFI mitigation. Put your own algorithm in. See below for instructions.


Use instructions:

 - psrenv preferred, or
 - /users/esmith/.conda/envs/py365 conda environment on green bank machine
 - type ' -h' to see help message

Inputs
------------
  -h, --help            show this help message and exit
  -i INFILE             String. Required. Name of input filename.
                        Automatically pulls from standard data directory. If
                        leading "/" given, pulls from given directory
  -rfi {SKurtosis,SEntropy,IQRM}
                        String. Required. RFI detection method desired.
  -m SK_M               Integer. Required. "M" in the SK equation. Number of
                        data points to perform SK on at once/average together
                        for spectrogram. ex. 1032704 (length of each block)
                        has prime divisors (2**9) and 2017. Default 512.
  -r {zeros,previousgood,stats}
                        String. Required. Replacement method of flagged data
                        in output raw data file. Can be
                        "zeros","previousgood", or "stats"
  -s SIGMA              Float. Sigma thresholding value. Default of 3.0 gives
                        probability of false alarm 0.001349
  -n N                  Integer. Number of inside accumulations, "N" in the SK
                        equation. Default 1.
  -v VEGAS_DIR          If inputting a VEGAS spectral line mode file, enter
                        AGBT19B_335 session number (1/2) and bank (C/D) ex
                        "1D".
  -newfile OUTPUT_BOOL  Copy the original data and output a replaced datafile.
                        Default True. Change to False to not write out a whole
                        new GUPPI file
  -d D                  Float. Shape parameter d. Default 1, but is different
                        in the case of low-bit quantization. Can be found (i
                        think) by running SK and changing d to be 1/x, where x
                        is the center of the SK value distribution.
  -npy RAWDATA          Boolean. True to save raw data to npy files. This is
                        storage intensive and unnecessary since blimpy.
                        Default is False
  -ms multiscale SK     String. Multiscale SK bin size. 
                        2 ints : Channel size / Time size, ex '-ms 42' Default '11'
  -mb mb		For loading multiple blocks at once. Helps with finding good
                        data for replacing flagged data, but can balloon RAM usage. 
                        Default 1.

#Assumes two polarizations
"""



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
from RFI_detection import *


#--------------------------------------
# Inputs
#--------------------------------------


#directories of interest
in_dir = '/data/rfimit/unmitigated/rawdata/'#leibniz only
out_dir = '/data/scratch/SKresults/'#leibniz only
jstor_dir = '/jetstor/scratch/SK_rawdata_results/'#leibniz only



parser = argparse.ArgumentParser(description="""function description""")


#input file
parser.add_argument('-i',dest='infile',type=str,required=True,help='String. Required. Name of input filename. Automatically pulls from standard data directory. If leading "/" given, pulls from given directory')

#replacement method
parser.add_argument('-r',dest='method',type=str,choices=['zeros','previousgood','stats','nans'], required=True,default='zeros',help='String. Required. Replacement method of flagged data in output raw data file. Can be "zeros","previousgood", nans or "stats"')

#write out a whole new raw file or just get SK/accumulated spectra results
parser.add_argument('-newfile',dest='output_bool',type=bool,default=True,help='Copy the original data and output a replaced datafile. Default True. Change to False to not write out a whole new GUPPI file')

#custom filename tag (for adding info not already covered in lines 187
parser.add_argument('-cust',dest='cust',type=str,default='',help='custom tag to add to end of filename')

#using multiple blocks at once to help stats replacement
parser.add_argument('-mult',dest='mb',type=int,default=1,help='load multiple blocks at once to help with stats/prevgood replacement')

#using multiple blocks at once to help stats replacement
parser.add_argument('-union',dest='union',type=int,default=1,help='Combine the polarizations in the flagging step. Default 1.')



#=================================================
# * * * * * * * * * * * * * *

#ID string of your RFI mitigation algorithm, for filenames and input parameters
IDstr = 'SK'

#parse the input arguments specific to your RFI mitigation algorithm

#example, for SK:


parser.add_argument('-SK_m',dest='SK_M',type=int,required=True,default=512,help='Integer. Required. "M" in the SK equation. Number of data points to perform SK on at once/average together for spectrogram. ex. 1032704 (length of each block) has prime divisors (2**9) and 2017. Default 512.')

#sigma thresholding
parser.add_argument('-SK_sigma',dest='sigma',type=float,default=3.0,help='Float. Sigma thresholding value. Default of 3.0 gives probability of false alarm 0.001349')

#number of inside accumulations, 'N' in the SK equation
parser.add_argument('-SK_n',dest='n',type=int,default=1,help='Integer. Number of inside accumulations, "N" in the SK equation. Default 1.')

#pick d in the case that it isn't 1. Required for low-bit quantization.
parser.add_argument('-SK_d',dest='d',type=float,default=1.,help='Float. Shape parameter d. Default 1, but is different in the case of low-bit quantization. Can be found (i think) by running SK and changing d to be 1/x, where x is the center of the SK value distribution.')

#multiscale bin shape.
parser.add_argument('-SK_ms',dest='ms',type=str,default='1,1',help='Multiscale SK. 2 ints : ChanSpec. Put a comma between. Default "1,1"')


args = parser.parse_args()
SK_M = args.SK_M
sigma = args.sigma
n = args.n
d = args.d
ms = (args.ms).split(',')
ms0 = int(ms[0])
ms1 = int(ms[1])

# * * * * * * * * * * * * * *
#=================================================

#load in the global arguments

infile = args.infile
method = args.method
#rawdata = args.rawdata
cust = args.cust
mb = args.mb
output_bool = args.output_bool
combine_flag_pols = args.union

#infile, method, rawdata, output_bool, cust, mb, combine_flag_pols = template_parse(parser)


#check infile, modify it to include in_dir if we don't give a full path to the file
infile = template_infile_mod(infile,in_dir)


#=================================================
# * * * * * * * * * * * * * *

#pattern for your parameters specific to the RFI
#example for SK:
outfile_pattern = "m{SK_M}_s{sigma}_ms{ms0}-{ms1}"


# any separate results filenames you need, in addition to the flags filename, put them here
npybase = out_dir+'npy_results/'+infile[len(in_dir):-4]


flags_filename = f"{npybase}_flags_{IDstr}_{outfile_pattern}_{cust}.npy"
spost_filename = f"{npybase}_spost_{IDstr}_{outfile_pattern}_{cust}.npy"


spect_filename = f"{npybase}_spect_{IDstr}_{outfile_pattern}_{cust}.npy"
sk_filename = f"{npybase}_skval_{IDstr}_{outfile_pattern}_{cust}.npy"
regen_filename = f"{npybase}_regen_{IDstr}_{outfile_pattern}_{cust}.npy"


#And then any one-off calculations at the beginning of the script

#threshold calc from sigma
SK_p = (1-scipy.special.erf(sigma/math.sqrt(2))) / 2
print('Probability of false alarm: {}'.format(SK_p))

#calculate thresholds
print('Calculating SK thresholds...')
lt, ut = SK_thresholds(SK_M, N = n, d = d, p = SK_p)
print('Upper Threshold: '+str(ut))
print('Lower Threshold: '+str(lt))

#calculate ms thresholds
ms_lt, ms_ut = SK_thresholds(SK_M*ms0*ms1, N = n, d = d, p = SK_p)
print('MS Upper Threshold: '+str(ms_ut))
print('MS Lower Threshold: '+str(ms_lt))



# * * * * * * * * * * * * * *
#=================================================


outfile = f"{jstor_dir}{infile[len(in_dir):-4]}_{IDstr}_{outfile_pattern}_mb{mb}_{cust}{infile[-4:]}"

if rawdata:
	print('Saving raw data to npy block style files')


#--------------------------------------
# Fun
#--------------------------------------


start_time = time.time()


if output_bool:
	template_check_outfile(infile,outfile)
	out_rawFile = open(outfile,'rb+')

#load file and copy
print('Opening file: '+infile)
rawFile = GuppiRaw(infile)


template_check_nblocks(rawFile,mb)


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
	print('Data shape: {} || block size: {}'.format(data.shape,data.nbytes))

	#save raw data?
	if rawdata:
		template_save_npy(data,block,npy_base)


#=================================================
# * * * * * * * * * * * * * *



	mismatch = num_timesamples % SK_M
	if (mismatch != 0):
		print('Warning: SK_M does not divide the amount of time samples')
		sys.exit()
	num_SKbins = int(num_timesamples/SK_M)

	#======================================
	#insert RFI detection of choice
	#
	#you currently have:
	#
	# data : 3D array of data from a single GUPPI/VPM block 'block'
	#	index [Channel,Spectrum,Polarization]
	#	these are np.complex64 complex channelized voltages
	#
	#	As well as any input parameters you specified above
	#
	#you need to output from this section:
	#
	# flag_chunk : 3D array of flags that correspond to the flags for this block
	#	shape can be the same as data or scrunched in the time axis
	#	(tscrunched flags are fast, so do that if you can)
	#	values follow the pattern
	#	0: unflagged || 1: flagged
	#
	# any intermediate output arrays 
	#	detection metrics, averaged spectra, etc. Make sure
	#	to fill them block-by-block using the if/else below
	#	and write them to disk at the end of the script
	#
	#
	#
	#	Ideally, your RFI code goes in RFI_detection_routines.py
	#	 and you simply place the function calls here.
	#
	#======================================

	


	#if you are making any intermediate numpy arrays (in addition to the flagging array), fill them here:

	if (block==0):

		flags_all = flag_chunk
	else:

		flags_all = np.concatenate((flags_all,flag_chunk),axis=1)

	# these will be written to disk at the end of the script

# * * * * * * * * * * * * * *
#=================================================

	#record flagging % in both polarizations
	flagged_pts_p1 += (1./numblocks) * ((100.*np.count_nonzero(flag_chunk[:,:,0]))/flag_chunk[:,:,0].size)
	flagged_pts_p2 += (1./numblocks) * ((100.*np.count_nonzero(flag_chunk[:,:,1]))/flag_chunk[:,:,1].size)


	#now flag shape is (chan,spectra,pol)
	#apply union of flags between the pols
	if combine_flag_pols:
		flag_chunk[:,:,0][flag_chunk[:,:,1]==1]=1
		flag_chunk[:,:,1][flag_chunk[:,:,0]==1]=1

	ts_factor = data.shape[1] // repl_chunk.shape[1]
	if (data.shape[1] % flag_chunk.shape[1] != 0):
		print('Flag chunk size is incompatible with block size')
		sys.exit()
		

	if method == 'zeros':
		#replace data with zeros
		data = repl_zeros(data,flag_chunk)

	if method == 'previousgood':
		#replace data with previous (or next) good
		data = previous_good(data,flag_chunk,ts_factor)

	if method == 'stats':
		#replace data with statistical noise derived from good datapoints
		data = statistical_noise_fir(data,flag_chunk,ts_factor)

	spost = template_averager(data,512)
	if (block==0):
		spost_all = spost
	else:
		spost_all = np.concatenate((spost_all,spost),axis=1)


	#Write back to copied raw file
	if output_bool:
		print('Re-formatting data and writing back to file...')
		for mb_i in range(mb):
			out_rawFile.seek(headersize,1)
			d1 = template_guppi_format(data[:,d1s*mb_i:d1s*(mb_i+1),:])
			out_rawFile.write(d1.tostring())




#save flags results
np.save(flags_filename,flags_all)
print(f'{flags_all.shape} Flags file saved to {flags_filename}')


#=================================================
# * * * * * * * * * * * * * *


#Any intermediate numpy arrays can be written out here, in addition to the flags array above




# * * * * * * * * * * * * * *
#=================================================


#tally up flags

template_print_flagstats(flags_all)


#clean up and end

print('Saved replaced data to '+outfile)


end_time = time.time()
elapsed = float(end_time-start_time)/60 

print('Program took {} minutes'.format(elapsed))

print('Done!')







