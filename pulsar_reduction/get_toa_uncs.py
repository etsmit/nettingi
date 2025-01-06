import os,sys
import numpy as np

basenm = "vegas_60299_75914_B0329+54_0003.0000"
basenm = "vegas_60299_75544_B0329+54_0001.0000"
basenm = "vegas_60299_76099_B0329+54_0004.0000"


unmit_toas = f"/data/scratch/SKresults/pulsar/{basenm}.raw/toas_new.zap.tim"
unmit_toas = f"/data/scratch/SKresults/pulsar/{basenm}.raw/toas_new.tim"
#unmit_toas = "/data/scratch/SKresults/pulsar/vegas_60299_75914_B0329+54_0003.0000.raw/toas_new.tim"
#unmit_toas = "/data/scratch/SKresults/pulsar/vegas_60299_75544_B0329+54_0001.0000.raw/toas_new.tim"
#unmit_toas = "/data/scratch/SKresults/pulsar/vegas_58966_25728_J1713+0747_0080.0000.raw/toas_new.tim"
#unmit_toas = "/data/scratch/SKresults/pulsar/vegas_58966_24913_J1713+0747_0076.0000.raw/toas_new.tim"
#unmit_toas = "/data/scratch/SKresults/pulsar/vegas_58966_24098_J1713+0747_0072.0000.raw/toas_new.tim"


mit_toas = "toas_new.tim"

def toa_reader(infile):
    f = open(infile,'r')
    lines = f.readlines()
    
    uncs = []
    mjds = []
    for i in range(0,len(lines),2):
        line = lines[i].split()
        mjds.append(float(line[2]))
        uncs.append(float(line[3]))
    mjds = np.array(mjds)
    uncs = np.array(uncs)
    return mjds,uncs

unmit_mjds,unmit_uncs = toa_reader(unmit_toas)
mit_mjds,mit_uncs = toa_reader(mit_toas)

unmit_unc = np.around( np.mean(unmit_uncs), 4 )
mit_unc = np.around( np.mean(mit_uncs), 4 )

print("Uncertainties")
print(f"Unmit: {unmit_unc} usec")
print(f"Mit:   {mit_unc} usec")
