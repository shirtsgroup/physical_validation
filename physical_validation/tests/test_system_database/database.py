###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
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
