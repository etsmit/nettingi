import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from argparse import ArgumentParser
from presto import sifting
import warnings
warnings.filterwarnings("ignore")

def read_accelcands(filenm):
    cands = []
    with open(filenm) as candfile:
        lines = candfile.readlines()
        line_to_read = 0
        while line_to_read < len(lines):
            line = lines[line_to_read]
            if line.startswith("#"):
                line_to_read += 1
            elif not line.startswith(" "):
                sline = line.strip().split()
                # candfiles have lines whose first column is formatted as
                # filename:candnum.  The join statement here ensures that
                # we capture the full filename if it happens to contain a ':'
                filename,candnum = (":".join(sline[0].split(":")[:-1]),
                                    int(sline[0].split(":")[-1]))
                DM = float(sline[1])
                SNR = float(sline[2])
                sigma = float(sline[3])
                numharm = int(sline[4])
                ipow = float(sline[5])
                cpow = float(sline[6])
                P = float(sline[7])*1e-3
                f = 1.0/P
                r = float(sline[8])
                T = r/f
                z = float(sline[9])
                # numhits is formatted as (numhits), so strip the parentheses
                numhits = int(sline[10].strip("()"))
                cand = sifting.Candidate(candnum,sigma,numharm,ipow,cpow,r,z,
                                         str(DM),filename,T)
                cand.snr = SNR
                
                if numhits > 1:
                    for nn in range(numhits):
                        line_to_read += 1
                        line = lines[line_to_read]
                        sline = line.split("=")
                        hit_DM = float(sline[1].split()[0])
                        hit_SNR = float(sline[2].split()[0])
                        hit_sigma = float(sline[3].split()[0])
                        cand.hits.append((hit_DM,hit_SNR,hit_sigma))
                cands.append(cand)
                line_to_read += 1

    return sifting.Candlist(cands)
                    

parser = ArgumentParser()
parser.add_argument("-m","--mitigated",required=True,
                    help="Mitigated .cands file")
parser.add_argument("-u","--unmitigated",required=True,
                    help="Unmitigated .cands file")
parser.add_argument("-p", "--period", required=True,type=float,
                    help="Pulsar period (seconds)")
parser.add_argument("-d", "--dm", required=True,type=float,
                    help="Pulsar DM (pc/cc)")
args = parser.parse_args()

center  = (1/args.period,args.dm)

mitigated_cands = read_accelcands(args.mitigated)
unmitigated_cands = read_accelcands(args.unmitigated)

with PdfPages(f"{os.path.basename(args.mitigated).replace('.cands','_cands_comparison.pdf')}") as pdf:
    fig1 = mitigated_cands.plot_summary()
    ax1 = fig1.gca()
    ax1.plot([1/args.period],[args.dm],color="C3",marker="o",mec="C3",
             mfc="none",ms=40,mew=5)
    ax1.set_title("Mitigated Candidates")
    plt.tight_layout()
    pdf.savefig(fig1)
    plt.show()

    fig2 = unmitigated_cands.plot_summary()
    ax2 = fig2.gca()
    ax2.plot([1/args.period],[args.dm],color="C3",marker="o",mec="C3",
             mfc="none",ms=40,mew=5)
    ax2.set_title("Unmitigated Candidates")
    plt.tight_layout()
    pdf.savefig(fig2)
    plt.show()

good_mitigated_cands = sorted(mitigated_cands.get_all_goodcands(),
                              key=lambda cand: cand.sigma,reverse=True)
good_unmitigated_cands = sorted(unmitigated_cands.get_all_goodcands(),
                                key=lambda cand: cand.sigma,reverse=True)

print(
f"""
+------------------+-------------------------+
|                  | Total Candidates        |
|------------------|-------------------------| 
|  Mitigated Data  | {len(mitigated_cands):^23.1f} |
|------------------|-------------------------|
| Unmitigated Data | {len(unmitigated_cands):^23.1f} |
+------------------+---------------------------+

"""
)

print(
f"""
+------------------+-------------------------+
|                  | Total 'Good' Candidates |
|------------------|-------------------------| 
|  Mitigated Data  | {len(good_mitigated_cands):^23.1f} |
|------------------|-------------------------|
| Unmitigated Data | {len(good_unmitigated_cands):^23.1f} |
+------------------+-------------------------+

"""
)

print("Top Ten Mitigated Candidates Sorted by Significance")
print("""+------------------+------------+---------+-----------+-----------+-------+
| Period (seconds) | DM (pc/cc) |  Sigma  |    SNR    |   Power   | Hits  |""")
for cand in good_mitigated_cands[:10]:
    print("|------------------|------------|---------|-----------|-----------|-------|")
    print((f"| {cand.p:^16.13f} | {cand.DM:^10.6f} | {cand.sigma:^7.2f} "
           f"| {cand.snr:^9.2f} | {cand.cpow:^9.2f} | {len(cand.hits):^5d} |"))
print("+------------------+------------+---------+-----------+-----------+-------+\n\n")

print("Top Ten Umitigated Candidates Sorted by Significance")
print("""+------------------+------------+---------+-----------+-----------+-------+
| Period (seconds) | DM (pc/cc) |  Sigma  |    SNR    |   Power   | Hits  |""")
for cand in good_unmitigated_cands[:10]:
    print("|------------------|------------|---------|-----------|-----------|-------|")
    print((f"| {cand.p:^16.13f} | {cand.DM:^10.6f} | {cand.sigma:^7.2f} "
           f"| {cand.snr:^9.2f} | {cand.cpow:^9.2f} | {len(cand.hits):^5d} |"))
print("+------------------+------------+---------+-----------+-----------+-------+")
