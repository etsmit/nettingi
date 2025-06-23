import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
parser.add_argument("-p","--period",required=True,
                    help="pulsar period")
parser.add_argument("-d","--dm",required=True,
                    help="dispersion measure")
parser.add_argument("-u","--unmit_dir",required=True,
                    help="unmitigated directory")   
args = parser.parse_args()

def bash(arg):
    print(arg)
    os.system(arg)

mit_mask_cands_fnm = f"{args.filename}_bary_search_mask.cands"
mit_nomask_cands_fnm = f"{args.filename}_bary_search_nomask.cands"

unmit_mask_cands_fnm = f"{args.unmit_dir}{args.filename}_bary_search_mask.cands"
unmit_nomask_cands_fnm = f"{args.unmit_dir}{args.filename}_bary_search_nomask.cands"

bash( f"python ACCEL_sift_mask.py > {mit_mask_cands_fnm}" )

bash( f"python ACCEL_sift_nomask.py > {mit_nomask_cands_fnm}" )

bash( "python cleancands.py" )

print(" --- No mask stats --- ")

bash( f"python compare_cands.py -u {unmit_nomask_cands_fnm} -m {mit_nomask_cands_fnm} -p {args.period} -d {args.dm}" )

print(" --- Mask stats --- ")

bash( f"python compare_cands.py -u {unmit_mask_cands_fnm} -m {mit_mask_cands_fnm} -p {args.period} -d {args.dm}" )







