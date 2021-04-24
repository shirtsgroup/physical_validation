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
Energy results of a system of 900 water molecules, used by the ensemble
check regression tests.
"""
import os

import numpy as np

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.UnitData.units("GROMACS"),
    system_data=pv.data.SystemData(
        natoms=900 * 3,
        nconstraints=900 * 3,
        ndof_reduction_tra=3,
        ndof_reduction_rot=0,
    ),
    ensemble={
        "NVT-low": pv.data.EnsembleData(
            ensemble="NVT",
            natoms=900 * 3,
            volume=3.01125 ** 3,
            temperature=298.15,
        ),
        "NVT-high": pv.data.EnsembleData(
            ensemble="NVT",
            natoms=900 * 3,
            volume=3.01125 ** 3,
            temperature=308.15,
        ),
        "NPT-lowT-lowP": pv.data.EnsembleData(
            ensemble="NPT", natoms=900 * 3, pressure=1.0, temperature=298.15
        ),
        "NPT-highT-lowP": pv.data.EnsembleData(
            ensemble="NPT", natoms=900 * 3, pressure=1.0, temperature=308.15
        ),
        "NPT-lowT-highP": pv.data.EnsembleData(
            ensemble="NPT", natoms=900 * 3, pressure=101.0, temperature=298.15
        ),
        "NPT-highT-highP": pv.data.EnsembleData(
            ensemble="NPT", natoms=900 * 3, pressure=101.0, temperature=308.15
        ),
    },
    description="900 water molecules, energy data only, flat and python data",
    simulation_keys="ensemble",
    observable_flat_file={
        "NVT-low": {
            "kinetic_energy": system_base_dir + "/NVT-low/kinetic.dat",
            "potential_energy": system_base_dir + "/NVT-low/potential.dat",
            "total_energy": system_base_dir + "/NVT-low/total.dat",
            "volume": system_base_dir + "/NVT-low/volume.dat",
        },
        "NVT-high": {
            "kinetic_energy": system_base_dir + "/NVT-high/kinetic.dat",
            "potential_energy": system_base_dir + "/NVT-high/potential.dat",
            "total_energy": system_base_dir + "/NVT-high/total.dat",
            "volume": system_base_dir + "/NVT-high/volume.dat",
        },
        "NPT-lowT-lowP": {
            "kinetic_energy": system_base_dir + "/NPT-lowT-lowP/kinetic.dat",
            "potential_energy": system_base_dir + "/NPT-lowT-lowP/potential.dat",
            "total_energy": system_base_dir + "/NPT-lowT-lowP/total.dat",
            "volume": system_base_dir + "/NPT-lowT-lowP/volume.dat",
        },
        "NPT-highT-lowP": {
            "kinetic_energy": system_base_dir + "/NPT-highT-lowP/kinetic.dat",
            "potential_energy": system_base_dir + "/NPT-highT-lowP/potential.dat",
            "total_energy": system_base_dir + "/NPT-highT-lowP/total.dat",
            "volume": system_base_dir + "/NPT-highT-lowP/volume.dat",
        },
        "NPT-lowT-highP": {
            "kinetic_energy": system_base_dir + "/NPT-lowT-highP/kinetic.dat",
            "potential_energy": system_base_dir + "/NPT-lowT-highP/potential.dat",
            "total_energy": system_base_dir + "/NPT-lowT-highP/total.dat",
            "volume": system_base_dir + "/NPT-lowT-highP/volume.dat",
        },
        "NPT-highT-highP": {
            "kinetic_energy": system_base_dir + "/NPT-highT-highP/kinetic.dat",
            "potential_energy": system_base_dir + "/NPT-highT-highP/potential.dat",
            "total_energy": system_base_dir + "/NPT-highT-highP/total.dat",
            "volume": system_base_dir + "/NPT-highT-highP/volume.dat",
        },
    },
    # This is the same data as from the flat files, but we're using numpy functionality
    # to read the data, so we can check that no error is introduced by our reading functions
    observable_as_array={
        "NVT-low": {
            "kinetic_energy": np.loadtxt(system_base_dir + "/NVT-low/kinetic.dat"),
            "potential_energy": np.loadtxt(system_base_dir + "/NVT-low/potential.dat"),
            "total_energy": np.loadtxt(system_base_dir + "/NVT-low/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NVT-low/volume.dat"),
        },
        "NVT-high": {
            "kinetic_energy": np.loadtxt(system_base_dir + "/NVT-high/kinetic.dat"),
            "potential_energy": np.loadtxt(system_base_dir + "/NVT-high/potential.dat"),
            "total_energy": np.loadtxt(system_base_dir + "/NVT-high/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NVT-high/volume.dat"),
        },
        "NPT-lowT-lowP": {
            "kinetic_energy": np.loadtxt(
                system_base_dir + "/NPT-lowT-lowP/kinetic.dat"
            ),
            "potential_energy": np.loadtxt(
                system_base_dir + "/NPT-lowT-lowP/potential.dat"
            ),
            "total_energy": np.loadtxt(system_base_dir + "/NPT-lowT-lowP/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NPT-lowT-lowP/volume.dat"),
        },
        "NPT-highT-lowP": {
            "kinetic_energy": np.loadtxt(
                system_base_dir + "/NPT-highT-lowP/kinetic.dat"
            ),
            "potential_energy": np.loadtxt(
                system_base_dir + "/NPT-highT-lowP/potential.dat"
            ),
            "total_energy": np.loadtxt(system_base_dir + "/NPT-highT-lowP/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NPT-highT-lowP/volume.dat"),
        },
        "NPT-lowT-highP": {
            "kinetic_energy": np.loadtxt(
                system_base_dir + "/NPT-lowT-highP/kinetic.dat"
            ),
            "potential_energy": np.loadtxt(
                system_base_dir + "/NPT-lowT-highP/potential.dat"
            ),
            "total_energy": np.loadtxt(system_base_dir + "/NPT-lowT-highP/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NPT-lowT-highP/volume.dat"),
        },
        "NPT-highT-highP": {
            "kinetic_energy": np.loadtxt(
                system_base_dir + "/NPT-highT-highP/kinetic.dat"
            ),
            "potential_energy": np.loadtxt(
                system_base_dir + "/NPT-highT-highP/potential.dat"
            ),
            "total_energy": np.loadtxt(system_base_dir + "/NPT-highT-highP/total.dat"),
            "volume": np.loadtxt(system_base_dir + "/NPT-highT-highP/volume.dat"),
        },
    },
)
