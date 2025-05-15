import numpy as np

import scipy as sp
import scipy.optimize
import scipy.special
import math as math
from blimpy import GuppiRaw

from .core import mitigateRFI

from .utils import *


class rfi_csp(mitigateRFI):
    #h
    def __init__(self, infile, repl_method, m, s, nbits=8, cust='', output_bool = True, mb=1, rawdata=False, ave_factor = 512):
        #user-given attributes
        self.det_method = 'CSP'
        self.repl_method = repl_method
        self.cust = cust
        self.output_bool = output_bool 
        self.mb = mb
        self.rawdata = rawdata
        self.ave_factor = ave_factor
        self.infile = infile 

        #default/hardcoded attributes

        