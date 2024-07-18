


import numpy as np


import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import pickle


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
        for i in range(len(infiles):
            freqs, data = load_tp(infiles[i])

    for i in range(len(spectra)):
        plt.plot(freqs, spectra[i], alpha=0.7, c=ten_clrs[i], label=labels=[i])

    plt.legend()
    plt.show()











