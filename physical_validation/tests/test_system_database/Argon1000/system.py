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
Trajectory of 1000 Argon atoms, used by the integrator convergence
regression tests.
"""
import os

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.UnitData.units("GROMACS"),
    system_data=pv.data.SystemData(
        natoms=1000,
        nconstraints=0,
        ndof_reduction_tra=3,
        ndof_reduction_rot=0,
    ),
    ensemble={
        "liquid": pv.data.EnsembleData(
            ensemble="NVE",
            natoms=1000,
            volume=3.60140 ** 3,
            energy=1,
        ),
    },
    description="1000 Argon atoms, different time step, flat file constant of motion only",
    simulation_keys="time step",
    time_step=[0.004, 0.002, 0.001, 0.0005, 0.00025, 0.000125],
    observable_flat_file={
        str(dt): {
            "constant of motion": f"{system_base_dir}/{dt}/constant_of_motion.dat"
        }
        for dt in [0.004, 0.002, 0.001, 0.0005, 0.00025, 0.000125]
    },
)
