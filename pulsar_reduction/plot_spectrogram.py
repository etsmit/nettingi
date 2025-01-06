import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob



sfile = glob.glob('*spect*.npy')[0]
ffile = glob.glob('*flags*.npy')[0]
mfile = glob.glob('*regen*.npy')[0]




s = np.load(sfile)
logs = np.log10(s)

f = np.load(ffile)
logsf = np.copy(logs)
logsf[f==1] = -3

sm = np.load(mfile)
logsm = np.log10(sm)


def implot(s,rf,bw,M,vmin=2.2,vmax=3.2):
	plt.figure(figsize=(14,10))
	ax=plt.axes()
	ax_lw = 3
	ax.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2)
	ax.spines['bottom'].set_linewidth(ax_lw)
	ax.spines['top'].set_linewidth(ax_lw)
	ax.spines['left'].set_linewidth(ax_lw)
	ax.spines['right'].set_linewidth(ax_lw)
	ms = s.shape[1]*(M*1e3*(s.shape[0]/(bw*1e6)))
	print(ms)
	#ext = [0,s.shape[1],1100,1900]
	ext = [0,ms/1e3,rf-0.5*bw,rf+0.5*bw]
	plt.imshow(s,interpolation='nearest',aspect='auto',cmap='hot',vmin=vmin,vmax=vmax,extent=ext)
	plt.ylabel('Frequency (MHz)',fontsize=20)
	plt.xlabel('Time (sec)',fontsize=20)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	cbar = plt.colorbar()
	for t in cbar.ax.get_yticklabels():
		t.set_fontsize(20)
	plt.tight_layout()
	plt.show()

print(s.shape)
print(np.count_nonzero(np.isnan(s)))


print(np.nanmean(s))
print(logs.shape)
print(np.count_nonzero(np.isnan(logs)))

print(np.nanmean(logs))

implot(logs[:,:,0], 1500, 800, 512)

implot(logsm[:,:,0], 1500, 800, 512)

implot(logsf[:,:,0], 1500, 800, 512)



