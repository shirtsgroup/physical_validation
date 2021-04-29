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
This file contains tests for the `physical_validation.data.trajectory` module.
"""
import numpy as np
import pytest

from ..data.trajectory_data import TrajectoryData
from ..util import error as pv_error


def test_trajectory_data_getters_and_setters() -> None:

    num_frames = 3
    num_atoms = 5
    position = np.random.random((num_frames, num_atoms, 3))
    velocity = np.random.random((num_frames, num_atoms, 3))

    # Check that objects created and populated in different ways are equivalent
    trajectory_data_1 = TrajectoryData()
    trajectory_data_1["position"] = position
    trajectory_data_1["velocity"] = velocity

    trajectory_data_2 = TrajectoryData()
    trajectory_data_2.velocity = velocity
    trajectory_data_2.position = position

    trajectory_data_3 = TrajectoryData(position=position, velocity=velocity)

    trajectory_data_2 = TrajectoryData()
    trajectory_data_2.velocity = velocity
    trajectory_data_2.position = position

    assert trajectory_data_1 == trajectory_data_2
    assert trajectory_data_1 == trajectory_data_3

    # Check setter error messages
    with pytest.raises(pv_error.InputError):
        # different number of frames
        trajectory_data_1.position = np.random.random((num_frames - 1, num_atoms, 3))
    with pytest.raises(pv_error.InputError):
        # different number of atoms
        trajectory_data_1.position = np.random.random((num_frames, num_atoms + 1, 3))
    with pytest.raises(pv_error.InputError):
        # different number of spatial dimensions
        trajectory_data_1.position = np.random.random((num_frames, num_atoms, 2))
    with pytest.raises(pv_error.InputError):
        # extra dimension in input
        trajectory_data_1.position = np.random.random((1, num_frames, num_atoms, 3))
    with pytest.raises(pv_error.InputError):
        # missing dimension in input
        trajectory_data_1.position = np.random.random(num_frames)
