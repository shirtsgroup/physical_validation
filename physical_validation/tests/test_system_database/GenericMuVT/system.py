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
Energy results of a generic muVT simulation, obtained by adding additional
independent species through bootstrapping of an actual simulation, and
using a very small chemical potential to avoid messing up the analysis.
Used to test the muVT ensemble check.
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
            ensemble="muVT", mu=[-0.0001, -37.5, -0.0001], volume=512, temperature=300
        ),
        "muVT-highT-lowMu": pv.data.EnsembleData(
            ensemble="muVT", mu=[-0.0001, -37.5, -0.0001], volume=512, temperature=305
        ),
        "muVT-lowT-highMu": pv.data.EnsembleData(
            ensemble="muVT", mu=[-0.0001, -37.0, -0.0001], volume=512, temperature=300
        ),
        "muVT-highT-highMu": pv.data.EnsembleData(
            ensemble="muVT", mu=[-0.0001, -37.0, -0.0001], volume=512, temperature=305
        ),
    },
    description="Generic muVT results, energy data only, flat files",
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
