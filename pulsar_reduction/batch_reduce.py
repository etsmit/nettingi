import os
import sys
import glob
from multiprocessing import Process, Queue, JoinableQueue

# dictionary of pulsar info
# Epoch/RA is key, value is list with [par file name, DM, period]
psr_info = {
    "B0329": ["B0329+54.par", 26.833, 0.71452],
    "B0355": ["B0355+54.par", 57.142, 0.156384],
    "J1713": ["J1713+0747.par", 15.992196, 0.0045701365],
}

pulsar = sys.argv[1]
basenm = sys.argv[2]
patt = sys.argv[3]


def bash(arg):
    try:
        os.sys(arg)
    except:  # noqa: E722
        print(f"Can't execute {arg}")


def put_q(q, a):
    q.put(bash)


my_q = JoinableQueue()
mp_list = []

dirs = glob.glob(patt)
for dir in dirs:
    print(dir)

input("Directories OK?")


for dir in dirs:
    parfile, dm, per = psr_info[pulsar]
    raw_filenm = ...
    arg = f"python reduce_raw_data.py -o {basenm} -d {dm} -p {parfile} "

    this_p = Process(target=bash, args=(arg,))
    mp_list.append(this_p)
    this_p.start()


for mp in mp_list:
    mp.join()
