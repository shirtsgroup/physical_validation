###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Pascal T. Merz <pascal.merz@me.com>                       #
#               Michael R. Shirts <michael.shirts@colorado.edu>           #
#                                                                         #
#    Copyright (c) 2017-2021 University of Colorado Boulder               #
#              (c) 2012      The University of Virginia                   #
#                                                                         #
###########################################################################
r"""
Physical validation suite for MD simulations

"""

__author__ = "Pascal T. Merz, and Michael R. Shirts"
__copyright__ = "2017-2021"
__credits__ = []
# TODO:
__license__ = "MIT"
__maintainer__ = "Pascal T. Merz, and Michael R. Shirts"
__email__ = "pascal.merz@me.com, michael.shirts@colorado.edu"

from . import data, ensemble, integrator, kinetic_energy, util
from ._version import get_versions

# Handle versioneer
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

__all__ = ["data", "ensemble", "integrator", "kinetic_energy", "util"]
