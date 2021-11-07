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
The system function allows to retrieve System objects from the database

When adding a new system to the database, this function must be extended
to encompass the new system.
"""
from typing import Dict

from .system import System


def system(system_name: str) -> System:
    if system_name == "Water900":
        from .Water900.system import system

        return system

    if system_name == "Octanol2":
        from .Octanol2.system import system

        return system

    if system_name == "Argon1000":
        from .Argon1000.system import system

        return system

    if system_name == "Water5":
        from .Water5.system import system

        return system

    if system_name == "Water300":
        from .Water300.system import system

        return system

    if system_name == "Octanol512":
        from .Octanol512.system import system

        return system

    if system_name == "DifluoromethaneGCMC":
        from .DifluoromethaneGCMC.system import system

        return system

    if system_name == "GenericMuVT":
        from .GenericMuVT.system import system

        return system

    raise KeyError(system_name)


def gromacs_files(system_name: str) -> Dict[str, Dict[str, str]]:
    if system_name == "Water5":
        from .Water5.system import gromacs_files

        return gromacs_files

    if system_name == "Octanol512":
        from .Octanol512.system import gromacs_files

        return gromacs_files

    raise KeyError(system_name)


def lammps_files(system_name: str) -> Dict[str, str]:
    if system_name == "Water300":
        from .Water300.system import lammps_files

        return lammps_files

    raise KeyError(system_name)
