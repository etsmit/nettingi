import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from presto import rfifind
import warnings
warnings.filterwarnings("ignore")


parser = ArgumentParser()
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
parser.add_argument("-u","--unmit_dir",required=True,
                    help="unmitigated directory")
args = parser.parse_args()

mitigated_fnm = f'{args.filename}_rfifind.mask'
unmitigated_fnm = f'{args.unmit_dir}{args.filename}_rfifind.mask'




mitigated_mask = rfifind.rfifind(mitigated_fnm)
mitigated_byte_mask = np.fromfile(
    mitigated_fnm.replace(".mask",".bytemask"),dtype=np.int8)
mitigated_byte_mask = mitigated_byte_mask.reshape(
    (mitigated_mask.nint,mitigated_mask.nchan)).astype(bool).astype(int)
unmitigated_mask = rfifind.rfifind(unmitigated_fnm)
unmitigated_byte_mask = np.fromfile(
    unmitigated_fnm.replace(".mask",".bytemask"),dtype=np.int8)
unmitigated_byte_mask = unmitigated_byte_mask.reshape(
    (unmitigated_mask.nint,unmitigated_mask.nchan)).astype(bool).astype(int)

if mitigated_mask.nint < unmitigated_mask.nint:
    mitigated_byte_mask = np.append(
        mitigated_byte_mask,np.zeros(
            (unmitigated_mask.nint-mitigated_mask.nint,mitigated_mask.nchan)),
        axis=0)
    print(("WARNING: Mitigated mask is "
           f"{unmitigated_mask.nint-mitigated_mask.nint} integrations shorter "
           "than unmitigated mask."))
elif mitigated_mask.nint > unmitigated_mask.nint:
    unmitigated_byte_mask = np.append(
        unmitigated_byte_mask,np.zeros(
            (mitigated_mask.nint-unmitigated_mask.nint,unmitigated_mask.nchan)),
        axis=0)
    print(("WARNING: Unmitigated mask is "
           f"{mitigated_mask.nint-unmitigated_mask.nint} integrations shorter "
           "than mitigated mask."))
else:
    pass

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.imshow(unmitigated_byte_mask-mitigated_byte_mask,origin="lower",
           aspect="auto",interpolation="hamming",cmap="bwr_r",vmin=-1,vmax=1)
ax1.set_xlabel("Channel")
ax1.set_ylabel("Index")
ax1.set_title("Blue is Good, Red is Bad")
ax2 = ax1.twinx().twiny()
ax2.plot(
    np.array([0,mitigated_mask.nchan-1])*mitigated_mask.df+mitigated_mask.lofreq,
    np.array([0,mitigated_mask.nint-1])*mitigated_mask.dtint,
    color="white",alpha=0)
ax2.set_xlabel("Freq (MHz)")
ax2.set_ylabel("Time (s)")
fig.suptitle("Mask Difference")
fig.tight_layout()
plt.savefig(
    f"{os.path.basename(mitigated_fnm).replace('.mask','_mask_comparison.png')}",
    dpi=300,bbox_inches="tight")
plt.show()
plt.close()

fig = plt.figure()
ax1a = fig.add_subplot(2,1,1)
ax1a.plot(100*mitigated_byte_mask.mean(axis=0),label="Mitigated",alpha=0.7)
ax1a.plot(100*unmitigated_byte_mask.mean(axis=0),label="Unmitigated",alpha=0.7)
ax1a.legend()
ax1a.set_xlabel("Channel")
ax1a.set_ylabel("Intervals Masked (%)")
ax1a.set_title("Masked Fraction Comparison")
ax1b = ax1a.twiny()
ax1b.plot(
    np.array([0,mitigated_mask.nchan-1])*mitigated_mask.df+mitigated_mask.lofreq,
    np.array([0,0]),color="white",alpha=1)
ax1b.set_xlabel("Freq (MHz)")
ax2a = fig.add_subplot(2,1,2)
ax2a.plot(100*mitigated_byte_mask.mean(axis=1))
ax2a.plot(100*unmitigated_byte_mask.mean(axis=1))
ax2a.set_xlabel("Interval")
ax2a.set_ylabel("Channels Masked (%)")
ax2b = ax2a.twiny()
ax2b.plot(
    np.array([0,mitigated_mask.nint-1])*mitigated_mask.dtint,np.array([0,0]),
          color="white",alpha=1)
ax2b.set_xlabel("Time (s)")
fig.tight_layout()
plt.savefig(
    f"{os.path.basename(mitigated_fnm).replace('.mask','_masked_fraction_comparison.png')}",
    dpi=300,bbox_inches="tight")
plt.show()
plt.close()

print(
f"""
+------------------+---------------------------+
|                  | Total Masked Fraction (%) |
|------------------|---------------------------| 
|  Mitigated Data  | {100*mitigated_byte_mask.mean():^25.1f} |
|------------------|---------------------------|
| Unmitigated Data | {100*unmitigated_byte_mask.mean():^25.1f} |
+------------------+---------------------------+
"""
)
