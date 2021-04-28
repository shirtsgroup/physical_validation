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
Trajectory of 512 united-atom octanol molecules, used by the kinetic
energy equipartition regression tests.
"""
import os

import numpy as np

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.UnitData.units("GROMACS"),
    system_data=pv.data.SystemData(
        natoms=512 * 10,
        nconstraints=512 * 9,
        ndof_reduction_tra=3,
        ndof_reduction_rot=0,
        mass=np.tile(
            [
                1.008,
                15.9994,
                14.027,
                14.027,
                14.027,
                14.027,
                14.027,
                14.027,
                14.027,
                15.035,
            ],
            512,
        ),
        molecule_idx=np.arange(0, 512 * 10, 10),
        nconstraints_per_molecule=9 * np.ones(512),
    ),
    ensemble={
        "GasPhase": pv.data.EnsembleData(
            ensemble="NVT",
            natoms=512 * 10,
            volume=1,
            temperature=298.15,
        ),
    },
    description="512 octanol molecules, system data and ensemble only",
    simulation_keys="ensemble",
)

gromacs_files = {
    "GasPhase": {
        "configuration": system_base_dir + "/gromacs_files/gas_sd.gro",
        "parameters": system_base_dir + "/gromacs_files/gas_sd.mdp",
        "topology": system_base_dir + "/gromacs_files/otl.top",
        "include_path": system_base_dir + "/gromacs_files/ff",
    }
}
