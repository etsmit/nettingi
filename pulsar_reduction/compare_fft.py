import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument("-m","--mitigated",required=True,
                    help="Mitigated FFT")
parser.add_argument("-u","--unmitigated",required=True,
                    help="Unmitigated FFT")
parser.add_argument("-p","--period",required=True,type=float,
                    help="Pulsar period (seconds")
parser.add_argument("-t", "--tsamp",required=True,type=float,
                    help="Sampling time (microseconds)")
args = parser.parse_args()

psr_freq = 1/args.period
mitigated_fft = np.fromfile(args.mitigated,dtype=np.complex64)
mitigated_freqs = np.fft.rfftfreq(
    2*len(mitigated_fft),d=args.tsamp*1e-6)[:-1]
unmitigated_fft = np.fromfile(args.unmitigated,dtype=np.complex64)
unmitigated_freqs = np.fft.rfftfreq(
    2*len(unmitigated_fft),d=args.tsamp*1e-6)[:-1]

plt.semilogx(mitigated_freqs,20*np.log10(np.abs(mitigated_fft)),alpha=0.7,
             label="Mitigated")
plt.semilogx(unmitigated_freqs,20*np.log10(np.abs(unmitigated_fft)),alpha=0.7,
             label="Unmitigated")
for n in range(16):
    label = "Pulsar w/ harmonics" if n == 0 else "_nolegend_"
    plt.axvline((n+1)*psr_freq,0.65,0.75,color="C3",label=label,zorder=-99)
plt.legend()
#plt.xlim(freqs.min(),freqs.max())
plt.xlabel("Freq (Hz)")
plt.ylabel("Power (dB)")
plt.savefig(
    f"{os.path.basename(args.mitigated).replace('.fft','_FFT_comparison.png')}",
    dpi=300,bbox_inches="tight")
plt.show()
plt.close()
