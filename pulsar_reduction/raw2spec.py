import os
import pickle
import numpy as np
from blimpy.guppi import GuppiRaw
#from pycycstat.utils import pfb
import scipy.signal as signal

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-r","--resolution",default=0.5,
                    help="Frequency resolution (kHz)")
parser.add_argument("-o","--outdir",default=".",
                    help="Output directory")
parser.add_argument("rawfiles",help="RAW files",nargs="*")
args = parser.parse_args()


def pfb(x, nchan, ntap, window="hann", fs=1.0, return_freqs=False, 
        force_complex=False):
    """
    Channelize data using a polyphase filterbank. From rlynch:pycycstat.

    Parameters
    ----------
    x : ndarray
       The input time series.
    nchan : int
       The number of channels to form.
    ntap : int
       The number of PFB taps to use.
    window : str
       The windowing function to use for the PFB coefficients.
    fs : float
       The sampling frequency of the input data.
    return_freqs : bool
       If True, return the center frequency of each channel.
    force_complex : bool
       If True, treat input as complex even if the imaginary component is zero.

    Returns
    -------
    x_pfb : ndarray
       The channelized data
    freqs : ndarray
       The center frequency of each channel.  Omitted if 
       return_freqs == False
    
    Notes
    -----
    If the input data are real-valued then only positive frequencies
    are returned.
    """
    real = np.isreal(x).all() and not force_complex
    h = signal.firwin(ntap*nchan,cutoff=1.0/nchan,window="rectangular")
    h *= signal.get_window(window,ntap*nchan)
    nwin = x.shape[0]//ntap//nchan
    x = x[:nwin*ntap*nchan].reshape((nwin*ntap,nchan)).T
    h = h.reshape((ntap,nchan)).T
    xs = np.zeros((nchan,ntap*(nwin-1)+1),dtype=x.dtype)
    for ii in range(ntap*(nwin-1)+1):
        xw = h*x[:,ii:ii+ntap]
        xs[:,ii] = xw.sum(axis=1)
    xs = xs.T
    xpfb = np.fft.rfft(xs,nchan,axis=1) if real else np.fft.fft(xs,nchan,axis=1)
    xpfb *= np.sqrt(nchan)

    if return_freqs:
        freqs = np.fft.rfftfreq(nchan,d=1.0/fs) if real else \
                np.fft.fftfreq(nchan,d=1.0/fs)
        return xpfb,freqs
    else:
        return xpfb


for rawfile in args.rawfiles:
    gr = GuppiRaw(rawfile)
    hdr0 = gr.read_first_header()
    fctr = hdr0["OBSFREQ"]
    bw = hdr0["OBSBW"]
    nchan = hdr0["OBSNCHAN"]
    chanbw = bw/nchan
    chanfreqs = fctr - 0.5*bw + chanbw*(np.arange(nchan)+0.5)
    nchan_pfb = 2**int(np.round(np.log2(np.abs(chanbw/(args.resolution/1e3)))))

    spectrum = np.zeros(nchan_pfb*nchan)

    gr.reset_index()
    for bb in range(gr.n_blocks):
        print("Working on block {0} of {1}".format(bb+1,gr.n_blocks))
        hdr,data = gr.read_next_data_block()
        x = data[:,:,0]
        y = data[:,:,1]
        for nn in range(data.shape[0]):
            xpfb = np.fft.fftshift(
                pfb(x[nn],nchan_pfb,12,force_complex=True),axes=-1)
            ypfb = np.fft.fftshift(
                pfb(y[nn],nchan_pfb,12,force_complex=True),axes=-1)
            spec = (np.mean(np.abs(xpfb)**2,axis=0)\
                    +np.mean(np.abs(ypfb)**2,axis=0))/2
            spectrum[nn*nchan_pfb:(nn+1)*nchan_pfb] += spec[::-1]
    spectrum /= gr.n_blocks
    pfb_chanbw = chanbw/nchan_pfb
    freqs = chanfreqs[0]-0.5*pfb_chanbw+np.arange(nchan_pfb*nchan)*pfb_chanbw
    bad_chans = np.arange(nchan)*nchan_pfb + nchan_pfb//2
    mask = np.zeros(spectrum.shape)
    mask[bad_chans] = 1
    masked_spectrum = np.ma.masked_array(spectrum,mask)
    out1 = (freqs,spectrum)
    out2 = (freqs,masked_spectrum)
    basenm = os.path.basename(rawfile)
    with open(os.path.join(args.outdir,basenm.replace(".raw",".spec.pkl")),"wb") as f: 
        pickle.dump(out1,f)
    with open(os.path.join(args.outdir,basenm.replace(".raw",".spec_mask.pkl")), "wb") as f: 
        pickle.dump(out2,f)



