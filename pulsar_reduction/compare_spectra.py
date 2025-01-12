import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
from argparse import ArgumentParser
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument("-m","--mitigated",required=True,
                    help="Mitigated spectrum")
parser.add_argument("-u","--unmitigated",required=True,
                    help="Unmitigated spectrum")
args = parser.parse_args()

with open(args.mitigated,"rb") as f: 
    mitigated_freqs,mitigated_spectrum = pickle.load(f)
with open(args.unmitigated,"rb") as f: 
    unmitigated_freqs,unmitigated_spectrum = pickle.load(f)

plt.plot(mitigated_freqs,mitigated_spectrum,alpha=0.7,label="Mitigated")
plt.plot(unmitigated_freqs,unmitigated_spectrum,alpha=0.7,label="Unmitigated")
plt.legend()
plt.xlabel("Freq (MHz)")
plt.ylabel("Counts")
plt.savefig(
    f"{os.path.basename(args.mitigated).replace('.spec_mask.pkl','_spectra_comparison.png')}",
    dpi=300,bbox_inches="tight")
plt.show()
plt.close()
