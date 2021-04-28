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
Energy and trajectory results of a system of 300 water molecules, used by the
LAMMPS parser tests.
"""
import os

import numpy as np

import physical_validation as pv

from ..system import System

system_base_dir = os.path.dirname(os.path.realpath(__file__))

system = System(
    units=pv.data.LammpsParser.units("real"),
    system_data=pv.data.SystemData(
        natoms=100 * 3,
        # LAMMPS parser doesn't set constraints
        nconstraints=0.0,
        ndof_reduction_tra=0,
        ndof_reduction_rot=0,
        mass=np.tile([15.9994, 1.008, 1.008], 100),
        molecule_idx=np.linspace(0, 100 * 3, 100, endpoint=False, dtype=int),
        # LAMMPS parser doesn't set constraints
        nconstraints_per_molecule=np.zeros(100),
    ),
    # Hack ensemble to also distinguish trajectories
    ensemble={
        "NVT": pv.data.EnsembleData(
            ensemble="NVT", natoms=100 * 3, volume=20 ** 3, temperature=300
        ),
    },
    description="100 water molecules, energy and trajectory data",
    simulation_keys="ensemble",
    time_step=[0],
    observable_flat_file={
        "NVT": {
            "kinetic_energy": system_base_dir + "/flat_files/kinetic.dat",
            "potential_energy": system_base_dir + "/flat_files/potential.dat",
            "total_energy": system_base_dir + "/flat_files/total.dat",
            "volume": system_base_dir + "/flat_files/volume.dat",
            "temperature": system_base_dir + "/flat_files/temperature.dat",
            "pressure": system_base_dir + "/flat_files/pressure.dat",
            "constant_of_motion": system_base_dir + "/flat_files/conserved.dat",
        }
    },
    trajectory_flat_file={
        "NVT": {
            "position": system_base_dir + "/flat_files/position.xyz",
            "velocity": system_base_dir + "/flat_files/velocity.xyz",
        }
    },
)

lammps_files = {
    "in": os.path.join(system_base_dir, "lammps_files", "water.in"),
    "log": os.path.join(system_base_dir, "lammps_files", "log.lammps"),
    "data": os.path.join(system_base_dir, "lammps_files", "water.lmp"),
    "dump": os.path.join(system_base_dir, "lammps_files", "dump.atom"),
}
