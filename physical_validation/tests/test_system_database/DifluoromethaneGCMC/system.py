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
Energy results of a grand canonical MC simulation of difluoromethane vapor,
used to test the muVT ensemble check.
"""
import os

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.UnitData.units("GROMACS"),
    system_data=pv.data.SystemData(),
    ensemble={
        "muVT-lowT-lowMu": pv.data.EnsembleData(
            ensemble="muVT", mu=-37.5, volume=512, temperature=300
        ),
        "muVT-highT-lowMu": pv.data.EnsembleData(
            ensemble="muVT", mu=-37.5, volume=512, temperature=305
        ),
        "muVT-lowT-highMu": pv.data.EnsembleData(
            ensemble="muVT", mu=-37.0, volume=512, temperature=300
        ),
        "muVT-highT-highMu": pv.data.EnsembleData(
            ensemble="muVT", mu=-37.0, volume=512, temperature=305
        ),
    },
    description="GCMC of difluoromethane molecules, energy data only, flat files",
    simulation_keys="ensemble",
    observable_flat_file={
        "muVT-lowT-lowMu": {
            "number_of_species": system_base_dir
            + "/muVT-lowT-lowMu/number_of_species.dat",
            "potential_energy": system_base_dir + "/muVT-lowT-lowMu/potential.dat",
        },
        "muVT-highT-lowMu": {
            "number_of_species": system_base_dir
            + "/muVT-highT-lowMu/number_of_species.dat",
            "potential_energy": system_base_dir + "/muVT-highT-lowMu/potential.dat",
        },
        "muVT-lowT-highMu": {
            "number_of_species": system_base_dir
            + "/muVT-lowT-highMu/number_of_species.dat",
            "potential_energy": system_base_dir + "/muVT-lowT-highMu/potential.dat",
        },
        "muVT-highT-highMu": {
            "number_of_species": system_base_dir
            + "/muVT-highT-highMu/number_of_species.dat",
            "potential_energy": system_base_dir + "/muVT-highT-highMu/potential.dat",
        },
    },
)
