


import numpy as np


import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import pickle
import glob  # noqa: F401


ten_clrs = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]


#infiles should be list of input pkl files. If position switched data, pair them up with 
#the ON scan first then the OFF scan.
#example for two position switched pairs, first one being the unmitigated:

#infiles = [
#    'vegas_60279_50072_Arp220_0004.0000.20.spec.pkl',
#    'vegas_60279_50334_Arp220_0005.0000.20.spec.pkl',
#    'vegas_60279_50072_Arp220_0004.0000_SK_m2032_ms1-1_thresh0.997_stats.20.spec_mask.pkl',
#    'vegas_60279_50334_Arp220_0005.0000_SK_m2032_ms1-1_thresh0.997_stats.20.spec_mask.pkl',
#    ]

#labels are just the list of custom labels you want to put in the legend
#labels = ['Unmitigated', 'Mitigated']

#and for this example, put ps=True so that it does the (on-off)/off

def pkl_plot(infiles, labels, ps=False):

    def pkll(inf):
        with open(inf,'rb') as f:                                         
            m_f,m_s = pickle.load(f)                                                  
        return m_f,m_s  

    def smooth(x,nchan):
        kernel = 1.0*np.ones(nchan) / nchan
        return np.convolve(x,kernel,mode='same')

    def load_tp(fname):
        freqs, data = pkll(fname)
        return freqs,data

    def load_ps_set(infiles):
        freqs, data_on = pkll(infiles[0])
        data_off = pkll(infiles[1])[1]
        psscan = (data_on-data_off)/data_off
        return freqs,psscan

    spectra = []

    if ps:
        for i in range(len(infiles)/2):
            freqs, data = load_ps_set(infiles[2*i],infiles[(2*i)+1])
            spectra.append(data)
            
    else:
        for i in range(len(infiles)):
            freqs, data = load_tp(infiles[i])
            spectra.append(data)

    for i in range(len(spectra)):
        plt.plot(freqs, spectra[i], alpha=0.7, c=ten_clrs[i], label=labels[i])

    plt.legend()
    plt.show()



#plot image of spectrogram
#s: 2D npy array you want to plot [chan, time]. NEEDS TO BE LOG10 SCALED
#    can be the original spectrogram, the mitigated spectrogram,
#    or you can combine the original and the flags (spect[flags==1]=0, s = spect)
#    to make the flagged portions black
#rf: MHz center frequency of the data (can read from the headers in the raw file)
#bw: MHz bandwidth of data (800 for the pulsar data)
#M:  averaging factor

#the latter 3 required arguments just help label the axes correctly
def implot(s,rf,bw,M,vmin=2.2,vmax=3.2):
    plt.figure(figsize=(14,10))
    ax=plt.axes()
    ax_lw = 3
    ax.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2)
    ax.spines['bottom'].set_linewidth(ax_lw)
    ax.spines['top'].set_linewidth(ax_lw)
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['right'].set_linewidth(ax_lw)
    ms = s.shape[1]*(M*1e3*(s.shape[0]/(bw*1e6)))
    print(ms)
    #ext = [0,s.shape[1],1100,1900]
    ext = [0,ms/1e3,rf-0.5*bw,rf+0.5*bw]
    plt.imshow(s,interpolation='nearest',aspect='auto',cmap='hot',vmin=vmin,vmax=vmax,extent=ext)
    plt.ylabel('Frequency (MHz)',fontsize=20)
    plt.xlabel('Time (sec)',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    cbar = plt.colorbar()
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(20)
    plt.tight_layout()
    plt.show()



def load_raw_flags(pattern,M):
    infiles = glob.glob(pattern)
    infiles.sort()

    for i,ff in enumerate(infiles):
        #print(ff, inflags[i])
        print(ff)
    print(f'{len(infiles)} files found')
    print('------------------------')

    input('Good?')


    init_f = np.load(infiles[0])
    out_f = np.empty((init_f.shape[0],init_f.shape[1]*len(infiles)/M,init_f.shape[2]))

    for i in range(len(infiles)):

        tf = np.load(infiles[i])
        a = np.reshape(tf,(tf.shape[0],-1,M,tf.shape[2]))

        nbins = tf.shape[1]//M
        start = i*nbins
        end = (i+1)*nbins

        out_f[:,start:end,:] = np.mean(a,axis=2)

    return out_f



















