#Support functions for mitigateRFI_template.py
#These should be used in all mitigateRFI variants

import numpy as np
import os
import sys

from numba import jit,prange
import nettingi
import glob

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
def template_bookkeeping(infile,out_patt,det):
    #outputs:
    # infile_raw_full  : full path and filename to input raw filename
    # output_raw_full  : full path and filename to output raw filename
    #                     makes output directory if needed
    # output_srdp_full : full path to SRDP files
    #                     makes the directory if needed
    #                     includes intermediate numpy files

#     #== input stuff ==
#     #these paths should all exist, so no mkdirs needed
#     #if they don't exist, 
    if infile[0] == "/":
        input_raw_dir_base = infile[:infile.rfind('/')+1]
    else:
        input_raw_dir_base = '/jetstor/scratch/rfimit/unmitigated/rawdata'
    infile_base = infile[infile.rfind('/')+1:infile.find('.')]
    print(infile,infile_base)
    
    #find input directory and file
    in_dir = f'{input_raw_dir_base}/{infile_base}/'
    infile_raw_full = in_dir+infile

#     #== output stuff ==

    #get output raw data file
    output_raw_dir_base = '/jetstor/scratch/rfimit/mitigated/rawdata'
    output_base = f'{infile_base}_{det}_{out_patt}'
    output_raw_dir = f'{output_raw_dir_base}/{output_base}/'
    if not os.path.exists(output_raw_dir):
        os.system(f'mkdir {output_raw_dir}')
    outfile_raw_full = f'{output_raw_dir}{infile[:-4]}_{det}_{out_patt}.raw'

    #get srdp results directory for mitigated file
    output_mit_srdp_dir_base = '/jetstor/scratch/rfimit/mitigated/reduced'
    if not os.path.exists(f'{output_mit_srdp_dir_base}/{output_base}'):
        os.system(f'mkdir {output_mit_srdp_dir_base}/{output_base}')
    output_mit_srdp_dir = f'{output_mit_srdp_dir_base}/{output_base}/{infile[:-4]}_{det}_{out_patt}/'
    if not os.path.exists(output_mit_srdp_dir):
        os.system(f'mkdir {output_mit_srdp_dir}')

    template_set_up_pulsar_reduction(output_mit_srdp_dir,infile_base)

    #get srdp results directory for unmitigated file
    # output_unmit_srdp_dir_base = '/jetstor/scratch/rfimit/unmitigated/reduced'
    # if not os.path.exists(f'{output_unmit_srdp_dir_base}/{output_base}'):
    #     os.system(f'mkdir {output_unmit_srdp_dir_base}/{output_base}')
    # output_unmit_srdp_dir = f'{output_unmit_srdp_dir_base}/{output_base}/{infile[:-4]}_{det}_{out_patt}/'
    # if not os.path.exists(output_unmit_srdp_dir):
    #     os.system(f'mkdir {output_unmit_srdp_dir}')

    print(infile_raw_full)
    print(outfile_raw_full)
    print(output_mit_srdp_dir)
    #print(output_unmit_srdp_dir)
    return infile_raw_full, outfile_raw_full, output_mit_srdp_dir

def template_set_up_pulsar_reduction(output_mit_srdp_dir, infile_base):
    psr_redux_dir = f"{nettingi.__file__[:-24]}pulsar_reduction/"
    psr_redux_scripts = glob.glob(f"{psr_redux_dir}*")
    for script in psr_redux_scripts:
        os.system(f"ln -s {script} {output_mit_srdp_dir}")
    psr = infile_base[18:23]
    if psr == 'J1713':
        os.system(f"ln -s {psr_redux_dir}parfiles/J1713+0747.par {output_mit_srdp_dir}")
    if psr == 'B0329':
        os.system(f"ln -s {psr_redux_dir}parfiles/B0329+54.par {output_mit_srdp_dir}")
    if psr == 'B0355':
        os.system(f"ln -s {psr_redux_dir}parfiles/B0355+54.par {output_mit_srdp_dir}")
    

#check that the outfile doesn't already exist, ask for overwrite confirmation 
def template_check_outfile(infile,outfile):
    print('Saving replaced data to '+outfile)
    print(infile,outfile)
    if os.path.isfile(outfile):
        #yn = input((f"The output file {outfile} already exists. Press 'y' to redo the copy, 'n' to continue without copying, or ctrl-c to end the script: "))
        yn='n'
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
    save_fname = npy_base+'_block'+block_fname+'.npy'
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
    #print(f'loading {hfile} for FIR coefficients')
    h = np.load(hfile)
    dec = h[::2*f.shape[0]]
    if ts_factor!=1:
        pulse = np.ones((1,ts_factor,1))
        f = np.kron(f,pulse)
    for pol in prange(f.shape[2]):
        for i in prange(f.shape[0]):
            #if ts_factor!=1:
            #    for tb in prange(f.shape[1]):
            #        if f[i,tb,pol] == 1:

                        #SK_M = ts_factor
                        #pulse = np.ones(ts_factor)
                        
 
            #            std_real,std_imag = adj_chan_good_data(a[:,tb*ts_factor:(tb+1)*ts_factor,pol],f[:,tb,pol],i)
                        
            #            (a[i,tb*SK_M:(tb+1)*SK_M,pol].real) = noise_filter(0,std_real,SK_M,dec)
                    
            #            (a[i,tb*SK_M:(tb+1)*SK_M,pol].imag) = noise_filter(0,std_imag,SK_M,dec)


            #else:
            bad_data_size = np.count_nonzero(f[i,:,pol])
            if bad_data_size > 0:
                #print(bad_data_size)
                std_real,std_imag = adj_chan_good_data(a[:,:,pol],f[:,:,pol],i)
                #print( noise_filter(0,std_real,bad_data_size,dec).shape)
                #print( a[i,:,pol].real[f[i,:,pol] == 1].shape)
                #print('--------------------')
                a[i,:,pol].real[f[i,:,pol] == 1] = noise_filter(0,std_real,bad_data_size,dec)[:bad_data_size]
                a[i,:,pol].imag[f[i,:,pol] == 1] = noise_filter(0,std_imag,bad_data_size,dec)[:bad_data_size]
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
    #print(a.shape,f.shape)
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



def statistical_noise_alt_fir(a,f,SK_M):
	"""
	Replace flagged data with statistical noise.
	- alt version that only takes the same time bin at adjacent-ish freq channels
	- fir version that filters noise by pfb coefficients to get desired spectral response
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
	print('stats....')
	print('fshape',f.shape)
	#find correct PFB coefficents
	#nchan = str(f.shape[0]*4).zfill(4)
	nchan = str(f.shape[0]).zfill(4)
	#print(nchan,type(nchan))	
	hfile = '/users/esmith/RFI_MIT/PFBcoeffs/c0800x'+nchan+'_x14_7_24t_095binw_get_pfb_coeffs_h.npy'
	h = np.load(hfile)
	dec = h[::2*f.shape[0]]

	for pol in prange(f.shape[2]):
		for i in prange(f.shape[0]):
			for tb in prange(f.shape[1]//SK_M):
				if f[i,SK_M*tb,pol] == 1:
			#find clean data points from same channel and polarization
			#good_data = a[i,:,pol][f[i,:,pol] == 0]
 
			#print(a[:,:,pol].shape,f[:,:,pol].shape)
					ave_real,ave_imag,std_real,std_imag = adj_chan_good_data_alt(a[:,tb*SK_M:(tb+1)*SK_M,pol],f[:,tb*SK_M:(tb+1)*SK_M,pol],i,SK_M,tb)
					#print('generating representative awgn..')
					(a[i,tb*SK_M:(tb+1)*SK_M,pol].real)[f[i,tb*SK_M:(tb+1)*SK_M,pol] == 1] = noise_filter(ave_real,std_real,SK_M,dec)
					
					(a[i,tb*SK_M:(tb+1)*SK_M,pol].imag)[f[i,tb*SK_M:(tb+1)*SK_M,pol] == 1] = noise_filter(ave_imag,std_imag,SK_M,dec)
	return a


def adj_chan_good_data_alt(a,f,c,SK_M,tb):
	"""
	Return mean/std derived from unflagged data in adjacent channels
	- alt version that only takes from adjacent channel at the SAME TIME BIN
	Parameters
	-----------
	a : ndarray
		3-dimensional array of original power values. Shape (Num Channels , Num Raw Spectra , Npol)
	f : ndarray
		3-dimensional array of flags. 1=RFI detected, 0 no RFI. Shape (Num Channels , Num Raw Spectra , Npol), should be same shape as a.
	c : int
		Channel of interest
	SK_M : int
		M parameter from SK equation
	tb : int
		time bin of interest
	
	Returns
	-----------
	ave_real : float
		average value of unflagged real data
	ave_imag : float
		average value of unflagged imaginary data
	std_real : float		
		standard deviation of unflagged real data
	std_imag : float
		standard deviation of unflagged imaginary data
	"""
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
	while (good_data.size==0):
		adj += 1
		if (c-adj >= 0):
			good_data = np.append(good_data,a[c-adj,:][f[c-adj,:] == 0])
		if (c+adj < a.shape[0]):
			good_data = np.append(good_data,a[c+adj,:][f[c+adj,:] == 0])
		#we gotta quit at some point, just give it flagged data from same channel
		#if we go 8% of the spectrum away
		if adj == int(a.shape[0]*0.08):
			good_data = a[c,:]
			break

	ave_real = np.mean(good_data.real)
	ave_imag = np.mean(good_data.imag)
	std_real = np.std(good_data.real)
	std_imag = np.std(good_data.imag)

	return ave_real, ave_imag, std_real,std_imag



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


def template_print_flagstats(flags_array,end):

    if end:
        add = '--Final--\n'
    else:
        add = ''
    print(f'{add}Pol 0: {np.around(100*np.mean(flags_array[:,:,0]),2)}% flagged')
    print(f'{add}Pol 1: {np.around(100*np.mean(flags_array[:,:,1]),2)}% flagged')

    uf = flags_array[:,:,0]
    uf[flags_array[:,:,1] == 1] = 1

    print(f'Union: {np.around(100*np.mean(uf),2)}% flagged')
    if end:
        print('--Final--')



#@jit(parallel=True)
def template_calc_ave(data,m):
    out = np.zeros((data.shape[0],data.shape[1]//m,data.shape[2]),dtype=np.float64)
    s = np.abs(data)**2
    #step1_p0 = np.ascontiguousarray(np.reshape(s[:,:,0], (s.shape[0],-1,m)))
    #step1_p1 = np.ascontiguousarray(np.reshape(s[:,:,1], (s.shape[0],-1,m)))
    #out[:,:,0] = np.nanmean(step1_p0,axis=2)
    #out[:,:,1] = np.nanmean(step1_p1,axis=2)
    a = np.reshape(s,(s.shape[0],-1,m,s.shape[2]))
    #numba nanmean cannot select by axis
    for pol in range(out.shape[2]):
        for chan in range(out.shape[0]):
            for tb in range(out.shape[1]):
                out[chan,tb,pol] = np.nanmean(a[chan,tb,:,pol])
    return out

def template_calc_std(data,m):
    out = np.zeros((data.shape[0],data.shape[1]//m,data.shape[2]))
    s = np.abs(data)**2

    step1_p0 = np.reshape(s[:,:,0], (s.shape[0],-1,m))
    step1_p1 = np.reshape(s[:,:,1], (s.shape[0],-1,m))
    out[:,:,0] = np.nanstd(step1_p0,axis=2)
    out[:,:,1] = np.nanstd(step1_p1,axis=2)
    return out

    






