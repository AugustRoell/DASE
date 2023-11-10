import sys
import platform
import os

from ._version import __version__

__all__ = [
    "config",
    "DASE",
    "materials",
    "profiles",
]


def version():
    r"""Version number of DASE.

    Returns:
      str: package version number
    """
    return __version__


def about():
    """Prints the installed version numbers for DASE, and some system info.

    Please include this information in bug reports.
    """

    print(
        "\nDASE is a software package for the simulation of radially symmetrical optical Schroedinger equation systems.",
    )

    print("Platform info:               {}".format(platform.platform()))
    print("Installation path:           {}".format(os.path.dirname(__file__)))
    print("Python version:              {}.{}.{}".format(*sys.version_info[0:3]))
    print("DASE version:                {}".format(__version__))
