import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from psrchive import Archive_load
from presto import fftfit
from presto.psr_utils import fft_rotate
from argparse import ArgumentParser
import warnings
warnings.filterwarnings("ignore")

def measure_phase(profile, template, rotate_prof=True):
    """
    measure_phase(profile, template):
        Call FFTFIT on the profile and template to determine the
            following parameters: shift,eshift,snr,esnr,b,errb,ngood
            (returned as a tuple).  These are defined as in Taylor's
            talk at the Royal Society.
    """
    c,amp,pha = fftfit.cprof(template)
    pha1 = pha[0]
    if (rotate_prof):
        pha = np.fmod(pha-np.arange(1,len(pha)+1)*pha1,2*np.pi)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit.fftfit(profile,amp,pha)
    return shift,eshift,snr,esnr,b,errb,ngood


parser = ArgumentParser()
#parser.add_argument("-m","--mitigated",required=True,
#                    help="Mitigated archive")
#parser.add_argument("-u","--unmitigated",required=True,
#                    help="Unmitigated archive")
parser.add_argument("-a", "--align",action="store_true",
                    help="Try to align profiles")
parser.add_argument("-f","--filename",required=True,
                    help="base filename")

args = parser.parse_args()

print(f'/jetstor/scratch/rfimit/unmitigated/reduced/{args.filename}.raw/{args.filename}_fold_sum.fits')

mitigated_archive = Archive_load(f'{args.filename}_fold_sum.fits')
unmitigated_archive = Archive_load(f'/jetstor/scratch/rfimit/unmitigated/reduced/{args.filename}.raw/{args.filename}_fold_sum.fits')

mitigated_archive.dedisperse()
unmitigated_archive.dedisperse()

mitigated_archive.remove_baseline()
unmitigated_archive.remove_baseline()

mitigated_archive.tscrunch()
unmitigated_archive.tscrunch()

mitigated_archive.fscrunch()
unmitigated_archive.fscrunch()

#mitigated_archive.bscrunch_to_nbin(512)
#unmitigated_archive.bscrunch_to_nbin(512)

mitigated_data = mitigated_archive.get_data().squeeze()
mitigated_data_tp = mitigated_data[0] + mitigated_data[1]
#mitigated_data_tp = np.roll(mitigated_data_tp,-1024)
#print(len(mitigated_data_tp))
unmitigated_data = unmitigated_archive.get_data().squeeze()
unmitigated_data_tp = unmitigated_data[0] + unmitigated_data[1]
#unmitigated_data_tp = np.roll(unmitigated_data_tp,-1024)
if args.align:
    shift = measure_phase(unmitigated_data_tp,mitigated_data_tp)[0]
    unmitigated_data_tp = fft_rotate(unmitigated_data_tp,-shift)
residual = mitigated_data_tp-unmitigated_data_tp
bins = np.arange(mitigated_archive.get_nbin())
phase = np.linspace(0,1,len(bins))

#calculate difference in peak flux / noise

#b0329
nx_start = 0
nx_end = int(0.8 * len(bins))

#j1713
#nx = np.r_[0:int(0.6 * len(bins)), int(0.9 * len(bins)):len(bins)]
#nx = np.r_[0:int(0.3 * len(bins)), int(0.6 * len(bins)):len(bins)]

#b0355
#nx = np.r_[0:int(0.6 * len(bins)), int(0.9 * len(bins)):len(bins)]
nx = np.r_[0:int(0.4 * len(bins)), int(0.6 * len(bins)):len(bins)]

#b0329
unmit_n = np.std(unmitigated_data_tp[nx_start:nx_end])
mit_n = np.std(mitigated_data_tp[nx_start:nx_end])

#j1713, b0355
#unmit_n = np.std(unmitigated_data_tp[nx])
#mit_n = np.std(mitigated_data_tp[nx])

print(unmit_n)
print(mit_n)

peak_diff = int(100.*np.around( np.max(mitigated_data_tp) / np.max(unmitigated_data_tp), 2))

mit_sn = np.max(mitigated_data_tp) / mit_n
unmit_sn = np.max(unmitigated_data_tp) / unmit_n

sn_diff = int(100.*np.around(mit_sn / unmit_sn,2))

print(f'Difference in peak power: {peak_diff} %')
print(f'Difference in S/N: {sn_diff} %')


c_red = '#FF0000'
c_bl = '#0000FF'
alph = 0.5
ax_lw=2
fontsz = 18

fig = plt.figure(figsize=(8,6))
ax = fig.gca()


ax.plot(phase,mitigated_data_tp,alpha=alph,label="Mitigated",c=c_red)
ax.plot(phase,unmitigated_data_tp,alpha=alph,label="Unmitigated",c=c_bl)
ax.axhline(0,c='k',linewidth=1)


ax.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2,labelsize=fontsz)
#ax.tick_params(axis='x',labelsize=1)
ax.spines['bottom'].set_linewidth(ax_lw)
ax.spines['top'].set_linewidth(ax_lw)
ax.spines['left'].set_linewidth(ax_lw)
ax.spines['right'].set_linewidth(ax_lw)






ax.legend(fontsize=fontsz)
ax.set_ylabel("Counts",fontsize=fontsz)
ax.set_xlabel("Pulse Phase",fontsize=fontsz)
plt.tight_layout()
plt.show()


#gs = plt.GridSpec(2,1,hspace=0,height_ratios=[4,1])
#fig = plt.figure()
#ax1 = fig.add_subplot(gs[0])
#ax1.plot(phase,mitigated_data_tp,alpha=0.5,label="Mitigated",c=c_red)
#ax1.plot(phase,unmitigated_data_tp,alpha=0.5,label="Unmitigated",c=c_bl)
#ax1.legend()
#ax1.set_ylabel("Counts")
#ax2 = fig.add_subplot(gs[1],sharex=ax1)
#ax2.plot(phase,residual,color="black")
#ax2.set_xlabel("Pulse Phase")
#ax2.set_ylabel("Residual")
#plt.savefig(f"{os.path.basename(args.mitigated).replace('.fits','_profile_comparison.png')}",dpi=300,bbox_inches="tight")
#plt.show()

