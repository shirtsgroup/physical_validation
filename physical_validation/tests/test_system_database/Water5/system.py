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
Energy and trajectory results of a system of 5 water molecules, used by the
GROMACS parser tests.
"""
import os

import numpy as np

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.UnitData.units("GROMACS"),
    system_data=pv.data.SystemData(
        natoms=5 * 3,
        nconstraints=5 * 3,
        ndof_reduction_tra=3,
        ndof_reduction_rot=0,
        mass=np.tile([15.9994, 1.008, 1.008], 5),
        molecule_idx=np.linspace(0, 5 * 3, 5, endpoint=False, dtype=int),
        nconstraints_per_molecule=3 * np.ones(5, dtype=int),
    ),
    # Ensemble is identical, but use this to distinguish trajectories
    ensemble={
        "full trajectory": pv.data.EnsembleData(
            ensemble="NPT", natoms=5 * 3, pressure=1, temperature=300
        ),
        "last frame only": pv.data.EnsembleData(
            ensemble="NPT", natoms=5 * 3, pressure=1, temperature=300
        ),
    },
    description="5 water molecules, energy and trajectory data",
    simulation_keys="ensemble",
    time_step=[0.001],
    observable_flat_file={
        "full trajectory": {
            "kinetic_energy": system_base_dir + "/flat_files/kinetic.dat",
            "potential_energy": system_base_dir + "/flat_files/potential.dat",
            "total_energy": system_base_dir + "/flat_files/total.dat",
            "volume": system_base_dir + "/flat_files/volume.dat",
            "temperature": system_base_dir + "/flat_files/temperature.dat",
            "pressure": system_base_dir + "/flat_files/pressure.dat",
            "constant_of_motion": system_base_dir + "/flat_files/conserved.dat",
        },
        "last frame only": {
            "kinetic_energy": system_base_dir + "/flat_files/kinetic.dat",
            "potential_energy": system_base_dir + "/flat_files/potential.dat",
            "total_energy": system_base_dir + "/flat_files/total.dat",
            "volume": system_base_dir + "/flat_files/volume.dat",
            "temperature": system_base_dir + "/flat_files/temperature.dat",
            "pressure": system_base_dir + "/flat_files/pressure.dat",
            "constant_of_motion": system_base_dir + "/flat_files/conserved.dat",
        },
    },
    trajectory_flat_file={
        "full trajectory": {
            "position": system_base_dir + "/flat_files/position_trajectory.xyz",
            "velocity": system_base_dir + "/flat_files/velocity_trajectory.xyz",
        },
        "last frame only": {
            "position": system_base_dir + "/flat_files/position_last_frame.xyz",
            "velocity": system_base_dir + "/flat_files/velocity_last_frame.xyz",
        },
    },
)

gromacs_files = {
    "parameters": system_base_dir + "/gromacs_files/mdout.mdp",
    "topology": system_base_dir + "/gromacs_files/topol.top",
    "final configuration": system_base_dir + "/gromacs_files/confout.gro",
    "energy": system_base_dir + "/gromacs_files/ener.edr",
    "trajectory": system_base_dir + "/gromacs_files/traj.trr",
}
