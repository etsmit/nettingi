import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from argparse import ArgumentParser
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument("-m","--mitigated",required=True,
                    help="Mitigated time series")
parser.add_argument("-u","--unmitigated",required=True,
                    help="Unmitigated time series")
args = parser.parse_args()

mitigated_time_series = np.fromfile(args.mitigated,dtype=np.float32)
unmitigated_time_series = np.fromfile(args.unmitigated,dtype=np.float32)

if len(mitigated_time_series) < len(unmitigated_time_series):
    diff = len(unmitigated_time_series) - len(mitigated_time_series)
    pad = np.ones(diff)*np.mean(mitigated_time_series)
    mitigated_time_series = np.append(mitigated_time_series,pad)
if len(mitigated_time_series) > len(unmitigated_time_series):
    diff = len(mitigated_time_series) - len(unmitigated_time_series)
    pad = np.ones(diff)*np.mean(unmitigated_time_series)
    unmitigated_time_series = np.append(unmitigated_time_series,pad)


fig = plt.figure()
gs = GridSpec(2,1,figure=fig,hspace=0.05,height_ratios=(2.,1))
ax1 = fig.add_subplot(gs[0])
ax1.plot(mitigated_time_series,alpha=0.7,label="Mitigated")
ax1.plot(unmitigated_time_series,alpha=0.7,label="Unmitigated")
ax1.set_ylabel("Power (arbitrary units)")
ax1.legend()
ax2 = fig.add_subplot(gs[1],sharex=ax1)
ax2.plot(
    (mitigated_time_series-unmitigated_time_series)/unmitigated_time_series*100)
ax2.set_xlabel("Sample")
ax2.set_ylabel("% Diff.")
ax1.set_title("Topocentric")
plt.savefig(f"{os.path.basename(args.mitigated).replace('.dat','_time_series_comparison.png')}",
            dpi=300,bbox_inches="tight")
plt.show()
plt.close()

print(
f"""
+------------------+---------------------------+
|                  | Time Series Standard Dev. |
|------------------|---------------------------| 
|  Mitigated Data  | {mitigated_time_series.std():^25.1f} |
|------------------|---------------------------|
| Unmitigated Data | {unmitigated_time_series.std():^25.1f} |
+------------------+---------------------------+
"""
)
