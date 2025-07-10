import matplotlib.pyplot as plt
import argparse


from presto.residuals import read_residuals


parser = argparse.ArgumentParser()
parser.add_argument("-f","--filename",required=True,
                    help="base filename")
# parser.add_argument("-u","--unmit_dir",required=True,
#                     help="unmitigated directory")
args = parser.parse_args()


unmit_red_dir = f'/jetstor/scratch/rfimit/unmitigated/reduced/{args.filename}.raw/'

unmit_resid = f"{unmit_red_dir}resid2.tmp"
mit_resid = "resid2.tmp"


c_red = '#FF0000'
c_bl = '#0000FF'
alph = 0.5
ax_lw=2
fontsz = 18

fig = plt.figure(figsize=(8,6))
ax = fig.gca()

ru = read_residuals(unmit_resid)

rm = read_residuals(mit_resid)

xu = ru.bary_TOA
yu = ru.postfit_phs
yu_err = ru.uncertainty

xm = rm.bary_TOA
ym = rm.postfit_phs
ym_err = rm.uncertainty


ax.errorbar(xu,yu,yerr=yu_err,marker='.',color=c_bl,linestyle='',label='Unmitigated',alpha=0.5,markersize=8)

ax.errorbar(xm,ym,yerr=ym_err,marker='.',color=c_red,linestyle='',label='Mitigated',alpha=0.5,markersize=8)
ax.legend()

ax.axhline(0,c='k')

ax.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2,labelsize=fontsz)
#ax.tick_params(axis='x',labelsize=1)
ax.spines['bottom'].set_linewidth(ax_lw)
ax.spines['top'].set_linewidth(ax_lw)
ax.spines['left'].set_linewidth(ax_lw)
ax.spines['right'].set_linewidth(ax_lw)


ax.legend(fontsize=fontsz)
ax.set_ylabel('Residuals (phase)',fontsize=fontsz)
ax.set_xlabel('DMJD',fontsize=fontsz)
plt.tight_layout()
plt.show()




