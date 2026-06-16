"""Top-level package for the project."""

from .core import *  # noqa: F403
from .sk import *  # noqa: F403
from .iqrmrfi import *  # noqa: F403
from .aof import *  # noqa: F403
from .mad import *  # noqa: F403
from .se import *  # noqa: F403
from .ConvRFIrfi import *  # noqa: F403
from .swnorm import *  # noqa: F403

__version__ = "0.1.0"

all = ["version"]


def version():
    """
    Version of the code

    :rtype: str
    """
    return __version__



