import os
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
parser.add_argument("-u","--unmit_dir",required=True,
                    help="unmitigated directory")
parser.add_argument("-p","--parfile",required=True,
                    help="parfile name")
args = parser.parse_args()

arg = f"paz -e zap.fits -r {args.filename}_fold_sum.fits"
print(arg)
os.system(arg)


arg = f"pam -e fscr --setnchn 8 {args.filename}_fold_sum.zap.fits"
print(arg)
os.system(arg)


arg1 = f"pat -F -f princeton -s {args.filename}_fold_sum.std {args.filename}_fold_sum.zap.fscr > toas.zap.tim"
#arg1 = "pat -F -f princeton -s vegas_58966_24098_J1713+0747_0072.0000_fold_sum.std vegas_58966_24098_J1713+0747_0072.0000_fold_sum.fscr > toas_new.tim"
print(arg1)
os.system(arg1)

#have to remove the last TOA, it is not calculated correctly in some cases.
f = open("toas.zap.tim",'r')
try:
    t = open('toas_new.zap.tim','x')
except:  # noqa: E722
    os.system('rm toas_new.zap.tim')
    t = open('toas_new.zap.tim','x')
lines = f.readlines()
for line in lines[:-1]:
    t.write(line)
f.close()
t.close()

arg2 = f"tempo -f {args.parfile} toas_new.zap.tim"
print(arg2)
os.system(arg2)

unmit_toas = f"{args.unmit_dir}toas_new.zap.tim"


mit_toas = "toas_new.zap.tim"

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
