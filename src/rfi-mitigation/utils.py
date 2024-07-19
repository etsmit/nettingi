#Support functions for mitigateRFI_template.py
#These should be used in all mitigateRFI variants

import numpy as np
import os,sys

import scipy as sp
import scipy.optimize
import scipy.special

import matplotlib.pyplot as plt

from numba import jit,prange

from tqdm import tqdm
import iqrm
import aoflagger

#from statsmodels.stats.weightstats import DescrStatsW



#get the template arguments
def template_parse(parser):

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

    #parse input variables
    args = parser.parse_args()
    infile = args.infile
    method = args.method
    rawdata = args.rawdata
    cust = args.cust
    mb = args.mb
    output_bool = args.output_bool
    combine_flag_pols = args.union
    return infile, method, rawdata, output_bool, cust, mb, combine_flag_pols


#modify input filename to include the right directory
def template_infile_mod(infile,in_dir):
    #input file
    #pulls from the raw data directory if full path not given
    if infile[0] != '/':
        infile = in_dir + infile
    else:
        in_dir = infile[:infile.rfind('/')+1]
        #infile = infile[infile.rfind('/')+1:]

    if infile[-4:] != '.raw':
        input("WARNING input filename doesn't end in '.raw'. Are you sure you want to use this file?")
    return infile


#set up all the input and output directories correctly
#set up all the input and output directories correctly
def template_bookkeeping(infile):

    #== input stuff ==
    #these paths should all exist, so no mkdirs needed
    #if they don't exist, 
    input_raw_dir_base = '/jetstor/scratch/rfimit/unmitigated/rawdata'
    infile_base = self.infile[:self.infile.find('.')]
    
    #find input directory and file
    in_dir = f'{input_raw_dir_base}/{infile_base}/'
    infile_raw_full = in_dir+self.infile

    #== output stuff ==
    output_raw_dir_base = '/jetstor/scratch/rfimit/mitigated/rawdata'
    output_base = f'{infile_base}_{self._output_pattern}'
    output_raw_dir = f'{output_raw_dir_base}/{output_base}/'
    if !os.path.exists(output_raw_dir):
        os.system(f'mkdir {output_raw_dir}')
    output_raw_full = f'{output_raw_dir}{infile[:-4]}_{self._output_pattern}.raw'


    output_srdp_dir_base = 'jetstor/scratch/rfimit/mitigated/reduced'
    output_srdp_dir = f'{output_srdp_dir_base}/{output_base}/{infile[:-4]}_{self._output_pattern}/'
    if !os.path.exists(output_srdp_dir):
        os.system(f'mkdir {output_srdp_dir}')
    return infile_raw_full, output_raw_full, output_srdp_full



#check that the outfile doesn't already exist, ask for overwrite confirmation 
def template_check_outfile(infile,outfile):
    print('Saving replaced data to '+outfile)
    print(infile,outfile)
    if os.path.isfile(outfile):
        yn = input((f"The output file {outfile} already exists. Press 'y' to redo the copy, 'n' to continue without copying, or ctrl-c to end the script: "))
        if yn=='y':
            print('Copying infile to outfile...')
            os.system('cp '+infile+' '+outfile)
    else:
        os.system('cp '+infile+' '+outfile)

#check that the number of blocks loaded at once is a divisible integer of the number of blocks in the file
def template_check_nblocks(rawFile,mb):
    numblocks = rawFile.find_n_data_blocks()
    print('File has '+str(numblocks)+' data blocks')
    #check for mismatched amount of blocks
    mismatch = numblocks % mb
    if (mismatch != 0):
        print(f'There are {numblocks} blocks and you set -mb {mb}, pick a divisible integer')
        sys.exit()

#read the first block's header
def template_print_header(rawFile):
    header,headersize = rawFile.read_header()
    print('Header size: {} bytes'.format(headersize))
    for line in header:
        print(line+':  '+str(header[line]))
    return headersize

#save numpy files
def template_save_npy(data,block,npy_base):
    block_fname = str(block).zfill(3)
    save_fname = npybase+'_block'+block_fname+'.npy'
    np.save(save_fname,data)






#=======================
#Replacement
#=======================

@jit
def repl_zeros(a,f):
    """
    Replace flagged data with 0's.

    Parameters
    -----------
    a : ndarray
        3-dimensional array of power values. Shape (Num Channels , Num Raw Spectra , Npol)
    f : ndarray
        3-dimensional array of flags. 1=RFI detected, 0 no RFI. Shape (Num Channels , Num Raw Spectra , Npol), should be same shape as a.
    
    
    Returns
    -----------
    out : ndarray
        3-dimensional array of power values with flagged data replaced. Shape (Num Channels , Num Raw Spectra , Npol)
    """
    ts = a.shape[1] // f.shape[1]
    if ts != 1:
        #for i in range(a.shape[1]):
        #    a[:,i,:][f[:,i//ts,:] == 1] = np.nan
        i = np.arange(a.shape[1])
        m = f[:,i//ts,:]
        a[m==1] = 0
    else:
        a[f==1] = 0
    return a



#@jit(nopython=False)
def repl_nans(a,f):
    """
    Replace flagged data with nans.

    Parameters
    -----------
    a : ndarray
        3-dimensional array of power values. Shape (Num Channels , Num Raw Spectra , Npol)
    f : ndarray
        3-dimensional array of flags. 1=RFI detected, 0 no RFI. Shape (Num Channels , Num Raw Spectra , Npol), should be same shape as a.
    
    
    Returns
    -----------
    out : ndarray
        3-dimensional array of power values with flagged data replaced. Shape (Num Channels , Num Raw Spectra , Npol)
    """
    ts = a.shape[1] // f.shape[1]
    if ts != 1:
        #for i in range(a.shape[1]):
        #    a[:,i,:][f[:,i//ts,:] == 1] = np.nan
        i = np.arange(a.shape[1])
        m = f[:,i//ts,:]
        a[m==1] = np.nan
    else:
        a[f==1]=np.nan
    return a

@jit
def repl_nans_jit(a,f):
    ts = a.shape[1] // f.shape[1]
    if ts != 1:
        i = np.arange(a.shape[1])
        m = f[:,i//ts,:]
        for pol in prange(a.shape[2]):
            for chan in prange(a.shape[0]):
                a[chan,:,pol][m[chan,:,pol]==1] = np.nan
    else:
        for pol in prange(a.shape[2]):
            for chan in prange(a.shape[0]):
                a[chan,:,pol][f[chan,:,pol]==1] = np.nan
    return a





#replace with statistical noise
# @jit(nopython=False)
def statistical_noise_fir(a,f,ts_factor):
    """
    Replace flagged data with statistical noise.
    - fir version that adds a fir in the noise
    Parameters
    -----------
    a : ndarray
        3-dimensional array of power values. Shape (Num Channels , Num Raw Spectra , Npol)
    f : ndarray
        3-dimensional array of flags. 1=RFI detected, 0 no RFI. Shape (Num Channels , Num Raw Spectra , Npol), should be same shape as a.

    
    Returns
    -----------
    out : np.random.normal(0,1,size=2048)ndarray
        3-dimensional array of power values with flagged data replaced. Shape (Num Channels , Num Raw Spectra , Npol)
    """
    #find correct PFB coefficents
    nchan = str(f.shape[0]).zfill(4)
    hfile = '/users/esmith/RFI_MIT/PFBcoeffs/c0800x'+nchan+'_x14_7_24t_095binw_get_pfb_coeffs_h.npy'
    print(f'loading {hfile} for FIR coefficients')
    h = np.load(hfile)
    dec = h[::2*f.shape[0]]
    if ts_factor!=1:
        pulse = np.ones((1,ts_factor,1))
        f = np.kron(f,pulse)
    for pol in prange(f.shape[2]):
        for i in prange(f.shape[0]):

                #for tb in prange(f.shape[1]):
                #    if f[i,tb,pol] == 1:

                #        SK_M = ts_factor
                #        pulse = np.ones(ts_factor)
                        
 
                #        std_real,std_imag = adj_chan_good_data(a[:,tb*SK_M:(tb+1)*SK_M,pol],f[:,tb,pol],i)
                        
                #        (a[i,tb*SK_M:(tb+1)*SK_M,pol].real) = noise_filter(0,std_real,SK_M,dec)
                    
                #        (a[i,tb*SK_M:(tb+1)*SK_M,pol].imag) = noise_filter(0,std_imag,SK_M,dec)


            #else:
            bad_data_size = np.count_nonzero(f[i,:,pol])
            if bad_data_size > 0:
                std_real,std_imag = adj_chan_good_data(a[:,:,pol],f[:,:,pol],i)

                a[i,:,pol].real[f[i,:,pol] == 1] = noise_filter(0,std_real,bad_data_size,dec)  
                a[i,:,pol].imag[f[i,:,pol] == 1] = noise_filter(0,std_imag,bad_data_size,dec)
    return a


# @jit
def adj_chan_good_data(a,f,c):
    """
    Return mean/std derived from unflagged data in adjacent channels 
    Parameters
    -----------
    a : ndarray
        3-dimensional array of original power values. Shape (Num Channels , Num Raw Spectra , Npol)
    f : ndarray
        3-dimensional array of flags. 1=RFI detected, 0 no RFI. Shape (Num Channels , Num Raw Spectra , Npol), should be same shape as a.
    c : int
        Channel of interest
    
    Returns
    -----------
    std_real : float        
        standard deviation of unflagged real data
    std_imag : float
        standard deviation of unflagged imaginary data
    """
    if len(f.shape) == 1:
        f = np.expand_dims(f,axis=1)
    #num_iter = num_iter
    #failed = failed
    #define adjacent channels and clear ones that don't exist (neg chans, too high)
    #adj_chans = [c-3,c-2,c-1,c,c+1,c+2,c+3]
    adj_chans = [c-1,c,c+1]
    adj_chans = [i for i in adj_chans if i>=0]
    adj_chans = [i for i in adj_chans if i<a.shape[0]]

    adj_chans = np.array(adj_chans,dtype=np.uint32)

    #set up array of unflagged data and populate it with any unflagged data from adj_chans channels
    good_data=np.empty(0,dtype=np.complex64)
    good_data = np.append(good_data,a[adj_chans,:][f[adj_chans,:] == 0])

    adj=1
    #keep looking for data in adjacent channels if good_data empty
    failed = 0    
    while (good_data.size==0):
        adj += 1
        if (c-adj >= 0):
            good_data = np.append(good_data,a[c-adj,:][f[c-adj,:] == 0])
        if (c+adj < a.shape[0]):
            good_data = np.append(good_data,a[c+adj,:][f[c+adj,:] == 0])
        #if we go 8% of the spectrum away, give up and give flagged data from same channel
        failed += 1
        if adj >= int(a.shape[0]*0.08):
            good_data = a[c,:]
            break

    std_real = np.std(good_data.real)
    std_imag = np.std(good_data.imag)

    return std_real,std_imag







@jit
def noise_filter(ave,std,msk,dec):
    """
    Create gaussian noise filtered by the correct PFB coefficients to mimic the VEGAS coarse channel SEFD
    Parameters
    -----------
    ave : float
        average/center value of intended noise 
    std : float
        standard deviation of noise (before FIR)
    msk : int
        M parameter of SK equation. Is also the amount of new data points to generate
    dec : decimated coefficient array to apply in FIR
    
    Returns
    -----------
    out_filtered : ndarray
        1-dimensional string of filtered gaussian noise to inject back over masked data
    """
    #make correctly scaled noise
    out = np.random.normal(ave,std,msk)
    #do FIR
    out_filtered = np.convolve(dec,out,mode='same')
    return out_filtered


@jit
def template_guppi_format(a):
    """
    takes array of np.complex64,ravels it and outputs as 1D array of signed 8 bit integers 
    ordered x1r,x1i,y1r,y1i,x2r,x2i,y2r,....
    Parameters
    -----------
    a : ndarray
        3-dimensional array of mitigated power values. Shape (Num Channels , Num Raw Spectra , Npol)
    Returns
    -----------
    out_arr : ndarray
        1-dimensional array of values to be written back to the copied data file
    """
    #init output
    out_arr = np.empty(shape=2*a.size,dtype=np.int8)
    #get real values, ravel, cast to int8
    arav = a.ravel()
    a_real = np.clip(np.floor(arav.real),-128,127).astype(np.int8)
    #get imag values, ravel, cast to int8
    a_imag = np.clip(np.floor(arav.imag),-128,127).astype(np.int8)
    #interleave
    out_arr[::2] = a_real
    out_arr[1::2] = a_imag
    return out_arr


def template_print_flagstats(flags_array):

    print(f'Pol 0: {np.around(100*np.mean(flags_array[:,:,0]),2)}% flagged')
    print(f'Pol 1: {np.around(100*np.mean(flags_array[:,:,1]),2)}% flagged')

    uf = flags_array[:,:,0]
    uf[flags_array[:,:,1] == 1] = 1

    print(f'Union: {np.around(100*np.mean(uf),2)}% flagged')



@jit(parallel=True)
def template_averager(data,m):
    out = np.zeros((data.shape[0],data.shape[1]//m,data.shape[2]),dtype=np.float64)
    s = np.abs(data)**2

    #step1_p0 = np.ascontiguousarray(np.reshape(s[:,:,0], (s.shape[0],-1,m)))
    #step1_p1 = np.ascontiguousarray(np.reshape(s[:,:,1], (s.shape[0],-1,m)))
    #out[:,:,0] = np.nanmean(step1_p0,axis=2)
    #out[:,:,1] = np.nanmean(step1_p1,axis=2)

    a = np.reshape(s,(s.shape[0],-1,m,s.shape[2]))
    #numba nanmean cannot select by axis
    for pol in prange(out.shape[2]):
        for chan in prange(out.shape[0]):
            for tb in prange(out.shape[1]):
                out[chan,tb,pol] = np.nanmean(a[chan,:,tb,pol])
    return out

def template_stdever(data,m):
    out = np.zeros((data.shape[0],data.shape[1]//m,data.shape[2]))
    s = np.abs(data)**2

    step1_p0 = np.reshape(s[:,:,0], (s.shape[0],-1,m))
    step1_p1 = np.reshape(s[:,:,1], (s.shape[0],-1,m))
    out[:,:,0] = np.nanstd(step1_p0,axis=2)
    out[:,:,1] = np.nanstd(step1_p1,axis=2)
    return out

def iqrm_power(data, radius, threshold):
    data = np.abs(data)**2
    flag_chunk = np.zeros(data.shape)
    for i in tqdm(range(data.shape[2])): # iterate through polarizations
        for j in range(data.shape[0]): # iterate through channels
            flag_chunk[j,:,i] = iqrm.iqrm_mask(data[j,:,i], radius = radius, threshold = threshold)[0]
    
#     avg_post = 
    return flag_chunk
    


def iqrm_std(data, radius, threshold, breakdown):
    """
    breakdown must be a factor of the time shape data[1].shape()
    """
# 	data_pol0 = stdever(np.abs(data[:,:,0])**2, breakdown) # make it a stdev
# # 	shape=np.expand_dims(shape, axis=2)
# 	flag_chunk = np.zeros((*data_pol0.shape[:2], 2))
# 	print('Data shape: {} || block size: {}'.format(flag_chunk.shape,flag_chunk.nbytes))
    
    data = template_stdever(np.abs(data)**2, breakdown)
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

    return flag_chunk

def iqrm_avg(data, radius, threshold, breakdown):
    """
    breakdown must be a factor of the time shape data[1].shape()
    """
    data = template_averager(np.abs(data)**2, breakdown)
    flag_chunk = np.zeros(data.shape)
    print('Flag shape: {} || block size: {}'.format(flag_chunk.shape,flag_chunk.nbytes))
    for i in tqdm(range(data.shape[2])): # iterate through polarizations
        for j in range(data.shape[0]): # iterate through channels
            flag_chunk[j,:,i] = iqrm.iqrm_mask(data[j,:,i], radius = radius, threshold = threshold)[0]

    return flag_chunk


def aof(data):
    nch = data.shape[0]
    ntimes = data.shape[1]
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
    tb_size = data.shape[1]//32    

    for tb in range(32):
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


