import os,sys
import numpy as np
import matplotlib.pyplot as plt

from presto.residuals import read_residuals

basenm = sys.argv[1]


unmit_resid = f"/jetstor/scratch/rfimit/unmitigated/reduced/{basenm}.raw/resid2.tmp"
mit_resid = "resid2.tmp"


c_red = '#FF0000'
c_bl = '#0000FF'

ru = read_residuals(unmit_resid)

rm = read_residuals(mit_resid)

xu = ru.bary_TOA
yu = ru.postfit_phs
yu_err = ru.uncertainty

xm = rm.bary_TOA
ym = rm.postfit_phs
ym_err = rm.uncertainty


plt.errorbar(xu,yu,yerr=yu_err,marker='.',color=c_bl,linestyle='',label='Unmitigated',alpha=0.5)

plt.errorbar(xm,ym,yerr=ym_err,marker='.',color=c_red,linestyle='',label='Mitigated',alpha=0.5)
plt.legend()

plt.axhline(0,c='k')

plt.ylabel('Residuals (phase)')
plt.xlabel('DMJD')
plt.show()




