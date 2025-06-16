import os,sys


import argparse
import time



parser = argparse.ArgumentParser()
parser.add_argument("-a", "--align",action="store_true",
                    help="Try to align profiles")
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
parser.add_argument("-p","--parfile",required=True,
                    help="parfile name")




args = parser.parse_args()
basenm = args.filename

unmit_red_dir = f'/jetstor/scratch/rfimit/unmitigated/reduced/{basenm}.raw/'

start_time = time.time

def bash(arg):
    print(arg)
    os.system(arg)

def waitfornextstep():
    input("Press Enter to continue...")


######################################
#  rfifind comparison
######################################


arg = f"python compare_rfi_masks.py -f {basenm}"
bash(arg)

waitfornextstep()

######################################
#  profile comparison
######################################

arg = f"python compare_profiles.py -f {basenm} -a"
bash(arg)

waitfornextstep()

######################################
#  TOAs comparison
######################################

arg = f"python get_TOAS.py {basenm} {}"
bash(arg)

waitfornextstep()



