import os
import argparse
import time
import numpy as np



parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
parser.add_argument("-p","--parfile",required=True,
                    help="parfile name")
parser.add_argument("-d","--dm",required=True,
                    help="dispersion measure")
parser.add_argument("-t","--period",required=True,
                    help="pulsar period")  



args = parser.parse_args()
basenm = args.filename

unmit_red_dir = f'/jetstor/scratch/rfimit/unmitigated/reduced/{basenm}.raw/'

def bash(arg):
    print(arg)
    os.system(arg)

def waitfornextstep():
    input("Press Enter to continue...")


######################################
#  rfifind comparison
######################################


arg = f"python compare_rfi_masks.py -f {basenm} -u {unmit_red_dir}"
bash(arg)

waitfornextstep()

######################################
#  profile comparison
######################################

arg = f"python compare_profiles.py -f {basenm} -u {unmit_red_dir} -a"
bash(arg)

waitfornextstep()

######################################
#  TOAs comparison
######################################

arg = f"python get_TOAs.py -f {basenm} -p {args.parfile} -u {unmit_red_dir}"
bash(arg)

waitfornextstep()

######################################
#  paz TOAs comparison
######################################

arg = f"python get_paz_TOAS.py -f {basenm} -p {args.parfile} -u {unmit_red_dir}"
bash(arg)

waitfornextstep()

######################################
#  candidate search comparison
######################################

#-p here is PERIOD, not PARFILE

arg = f"python do_cands.py -f {basenm} -p {args.period} -u {unmit_red_dir} -d {args.dm}"
bash(arg)



