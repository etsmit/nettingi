import sys,os



basenm = sys.argv[1]
per = sys.argv[2]
dm = sys.argv[3]


def bash(arg)
    print(arg)
    os.system(arg)



bash( f"python ACCEL_sift_mask.py > {basenm}_bary_search_mask.cands" )

bash( f"python ACCEL_sift_nomask.py > {basenm}_bary_search_nomask.cands" )

bash( "python cleancands.py" )

print(" --- No mask stats --- ")

bash( f"python compare_cands.py -u {unmit_nomask_cands} -m {mit_nomask_cands}" )









