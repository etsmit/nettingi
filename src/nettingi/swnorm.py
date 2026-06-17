# spectral entropy

import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
from blimpy import GuppiRaw

from .core import mitigateRFI

from .utils import template_bookkeeping


class rfi_swnorm(mitigateRFI):
    # h
    def __init__(
        self,
        infile,
        repl_method,
        m=512,
        alpha=1e-3,
        cust="",
        output_bool=True,
        mb=1,
        rawdata=False,
        ave_factor=512,
    ):
        # user-given attributes
        self.det_method = "SWNORM"
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile

        # default/hardcoded attributes

        self.m = m
        self.alpha = alpha

        self._outfile_pattern = f"m{self.m}_a{self.alpha}"

        self.infile_raw_full, self.outfile_raw_full, self.output_mit_srdp_dir = (
            template_bookkeeping(self.infile, self._outfile_pattern, self.det_method)
        )
        self._rawFile = GuppiRaw(self.infile_raw_full)
        # any separate results filenames you need, in addition to the flags filename, put them here
        self.npybase = self.infile[:-4]

        self._flags_filename = f"{self.output_mit_srdp_dir}{self.npybase}_flags_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._spect_filename = f"{self.output_mit_srdp_dir}{self.npybase}_spect_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._regen_filename = f"{self.output_mit_srdp_dir}{self.npybase}_regen_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        # self._outfile = f"{self._jetstor_dir}{infile[:-4]}_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_mb{self.mb}_{self.cust}{infile[-4:]}"

        # any derived thresholds/arrays
        self._ptest_filename = f"{self.output_mit_srdp_dir}{self.npybase}_ptest_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"
        self._stat_filename = f"{self.output_mit_srdp_dir}{self.npybase}_stat_{self.det_method}_{self.repl_method}_{self._outfile_pattern}_{self.cust}.npy"

        out = f"""input: {self.infile_raw_full}\noutput: {self.outfile_raw_full}\nspect: {self._spect_filename}"""
        print(out)

    def swnorm_detection(self, data):

        a = np.reshape(data, (data.shape[0], -1, self.m, data.shape[2]))

        # yippee, a scipy stats function for free, for me!
        print("swnorm...")
        shap = sp.stats.shapiro(a, axis=2, nan_policy="omit")
        ptest = shap.pvalue
        stat = shap.statistic

        # print(ptest.shape)

        print("flagging")
        # flag
        flags_block = np.zeros(ptest.shape, dtype=np.int8)
        flags_block[np.isnan(ptest)] = 1
        flags_block[ptest < self.alpha] = 1

        return flags_block, ptest, stat
