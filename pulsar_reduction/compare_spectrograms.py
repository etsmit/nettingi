import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from argparse import ArgumentParser
from presto.filterbank import FilterbankFile
import warnings
warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument("-m","--mitigated",required=True,nargs="*",
                    help="Mitigated filterbank files")
parser.add_argument("-u","--unmitigated",required=True,nargs="*",
                    help="Unmitigated filterbank files")
parser.add_argument("-s", "--sample",required=True,type=float,
                    help="Center sample number of range to extract")
parser.add_argument("-t", "--time", default=14.0, type=float,
                    help="Duration to plot (seconds)")
parser.add_argument("-c", "--channels", help="Channels to plot.  Range separated by '-', e.g. 128-132")
parser.add_argument("--smooth", action="store_true", help="Smooth the bandpass")
args = parser.parse_args()

mitigated_spectrogram = []
samples_read = samples_skipped = 0
for filterbank_file in args.mitigated:
    fb = FilterbankFile(filterbank_file)
    samples_to_extract = int(np.ceil(args.time/fb.dt/2)*2)
    start_sample = int(args.sample) - samples_skipped - samples_to_extract//2
    if start_sample < 0: start_sample = 0
    end_sample = int(args.sample) - samples_skipped + samples_to_extract//2
    if end_sample < fb.nspec:
        spectra = fb.get_spectra(start_sample,samples_to_extract-samples_read)
        for nn in range(spectra.numspectra):
            spec = spectra.get_spectrum(nn)[::-1]
            if args.smooth: spec -= medfilt(spec,11)
            mitigated_spectrogram.append(spec)
        break
    elif start_sample < fb.nspec and end_sample > fb.nspec:
        spectra = fb.get_spectra(start_sample,fb.nspec-start_sample)
        for nn in range(spectra.numspectra):
            spec = spectra.get_spectrum(nn)[::-1]
            if args.smooth: spec -= medfilt(spec,11)
            mitigated_spectrogram.append(spec)
        samples_read = int(fb.nspec - start_sample)
        samples_skipped += int(fb.nspec)
    else:
        samples_skipped += int(fb.nspec)
mitigated_spectrogram = np.array(mitigated_spectrogram)
freqs = fb.frequencies[::-1]

unmitigated_spectrogram = []
samples_read = samples_skipped = 0
total_samps = 0
for filterbank_file in args.unmitigated:
    fb = FilterbankFile(filterbank_file)
    total_samps += fb.nspec
    samples_to_extract = int(np.ceil(args.time/fb.dt/2)*2)
    start_sample = int(args.sample) - samples_skipped - samples_to_extract//2
    if start_sample < 0: start_sample = 0
    end_sample = int(args.sample) - samples_skipped + samples_to_extract//2
    if end_sample < fb.nspec:
        spectra = fb.get_spectra(start_sample,samples_to_extract-samples_read)
        for nn in range(spectra.numspectra):
            spec = spectra.get_spectrum(nn)[::-1]
            if args.smooth: spec -= medfilt(spec,11)
            unmitigated_spectrogram.append(spec)
        break
    elif start_sample < fb.nspec and end_sample > fb.nspec:
        spectra = fb.get_spectra(start_sample,fb.nspec-start_sample)
        for nn in range(spectra.numspectra):
            spec = spectra.get_spectrum(nn)[::-1]
            if args.smooth: spec -= medfilt(spec,11)
            unmitigated_spectrogram.append(spec)
        samples_read = int(fb.nspec - start_sample)
        samples_skipped += int(fb.nspec)
    else:
        samples_skipped += int(fb.nspec)
unmitigated_spectrogram = np.array(unmitigated_spectrogram)

if args.channels is not None:
    low_chan,high_chan = map(int, args.channels.split("-"))
    low_freq = freqs[low_chan]
    high_freq = freqs[high_chan]
    mitigated_spectrogram = mitigated_spectrogram[:,low_chan:high_chan]
    unmitigated_spectrogram = unmitigated_spectrogram[:,low_chan:high_chan]
else:
    low_freq = freqs[0]
    high_freq = freqs[-1]

low_samp = args.sample-args.time/fb.dt//2
if low_samp < 0: low_samp = 0
high_samp = args.sample + args.time/fb.dt//2
if high_samp > total_samps: high_samp = total_samps
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax1.imshow(mitigated_spectrogram,interpolation="none",aspect="auto",
           origin="lower",extent=[low_freq,high_freq,low_samp,high_samp])
ax1.set_xlabel("Frequency (MHz)")
ax1.set_ylabel("Sample Number")
ax1.set_title("Mitigated")

ax2 = fig.add_subplot(1,2,2,sharex=ax1,sharey=ax1)
ax2.imshow(unmitigated_spectrogram,interpolation="none",aspect="auto",
           origin="lower",extent=[low_freq,high_freq,low_samp,high_samp])
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_xlabel("Frequency (MHz)")
ax2.set_title("Unmitigated")
fig.savefig(f"{os.path.basename(args.mitigated[0]).replace('.fil','_spectrograms_comparison_sample{args.sample}.png')}",dpi=300,bbox_inches="tight")
plt.show()
plt.close()
