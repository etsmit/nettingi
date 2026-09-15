"""Top-level package for the project."""

from .aof import *
from .ConvRFIrfi import *
from .core import *
from .iqrmrfi import *
from .mad import *
from .se import *
from .sk import *
from .swnorm import *

__version__ = "0.1.0"

all = ["version"]


def version():
    """
    Version of the code

    :rtype: str
    """
    return __version__
