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
import itertools
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
    # Hack ensemble to also distinguish trajectories
    ensemble={
        "NVT full trajectory": pv.data.EnsembleData(
            ensemble="NVT", natoms=5 * 3, volume=1.86210 ** 3, temperature=300
        ),
        "NVT last frame only": pv.data.EnsembleData(
            ensemble="NVT", natoms=5 * 3, volume=1.86210 ** 3, temperature=300
        ),
        "NPT full trajectory": pv.data.EnsembleData(
            ensemble="NPT", natoms=5 * 3, pressure=1, temperature=300
        ),
        "NPT last frame only": pv.data.EnsembleData(
            ensemble="NPT", natoms=5 * 3, pressure=1, temperature=300
        ),
    },
    description="5 water molecules, energy and trajectory data",
    simulation_keys="ensemble",
    time_step=[0.001],
    observable_flat_file={
        f"{ensemble} {trajectory}": {
            "kinetic_energy": system_base_dir + f"/flat_files/{ensemble}/kinetic.dat",
            "potential_energy": system_base_dir
            + f"/flat_files/{ensemble}/potential.dat",
            "total_energy": system_base_dir + f"/flat_files/{ensemble}/total.dat",
            "volume": system_base_dir + f"/flat_files/{ensemble}/volume.dat",
            "temperature": system_base_dir + f"/flat_files/{ensemble}/temperature.dat",
            "pressure": system_base_dir + f"/flat_files/{ensemble}/pressure.dat",
            "constant_of_motion": system_base_dir
            + f"/flat_files/{ensemble}/conserved.dat",
            "number_of_species": system_base_dir
            + f"/flat_files/{ensemble}/number_of_species.dat",
        }
        for ensemble, trajectory in itertools.product(
            ["NPT", "NVT"], ["full trajectory", "last frame only"]
        )
    },
    trajectory_flat_file={
        **{
            f"{ensemble} full trajectory": {
                "position": system_base_dir
                + f"/flat_files/{ensemble}/position_trajectory.xyz",
                "velocity": system_base_dir
                + f"/flat_files/{ensemble}/velocity_trajectory.xyz",
            }
            for ensemble in ["NPT", "NVT"]
        },
        **{
            f"{ensemble} last frame only": {
                "position": system_base_dir
                + f"/flat_files/{ensemble}/position_last_frame.xyz",
                "velocity": system_base_dir
                + f"/flat_files/{ensemble}/velocity_last_frame.xyz",
            }
            for ensemble in ["NPT", "NVT"]
        },
    },
)

gromacs_files = {
    ensemble: {
        "configuration": system_base_dir + "/gromacs_files/conf.gro",
        "parameters": system_base_dir + f"/gromacs_files/{ensemble}/mdout.mdp",
        "topology": system_base_dir + "/gromacs_files/topol.top",
        "final configuration": system_base_dir
        + f"/gromacs_files/{ensemble}/confout.gro",
        "energy": system_base_dir + f"/gromacs_files/{ensemble}/ener.edr",
        "trajectory": system_base_dir + f"/gromacs_files/{ensemble}/traj.trr",
    }
    for ensemble in ["NPT", "NVT"]
}
