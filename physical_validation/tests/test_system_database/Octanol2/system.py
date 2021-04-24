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
Trajectory of two united-atom octanol molecules, used by the kinetic
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
        natoms=2 * 10,
        nconstraints=2 * 9,
        ndof_reduction_tra=3 / (512 / 2),
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
            2,
        ),
        molecule_idx=np.arange(0, 2 * 10, 10),
        nconstraints_per_molecule=9 * np.ones(2),
    ),
    ensemble={
        "GasPhase": pv.data.EnsembleData(
            ensemble="NVT",
            natoms=2 * 10,
            volume=1,
            temperature=298.15,
        ),
    },
    description="2 octanol molecules, flat file trajectories only",
    simulation_keys="ensemble",
    trajectory_flat_file={
        "GasPhase": {
            "position": system_base_dir + "/GasPhase/position.xyz",
            "velocity": system_base_dir + "/GasPhase/velocity.xyz",
        }
    },
)
