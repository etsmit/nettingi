#getcands.py

import sys,os

arg = sys.argv[1]

print('accel sift mask...')
os.system(f'python ACCEL_sift_mask.py > {arg}_bary_search_mask.cands')
print('accel sift nomask...')
os.system(f'python ACCEL_sift_nomask.py > {arg}_bary_search_nomask.cands')

print('cleaning...')
os.system('python cleancands.py')





