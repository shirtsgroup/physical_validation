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
This file contains tests for the `physical_validation.data.observable_data` module.
"""
import numpy as np
import pytest

from ..data.observable_data import ObservableData
from ..util import error as pv_error


def test_observable_data_getters_and_setters() -> None:

    # Check that newly created observable data object has `None` frames
    observable_data = ObservableData()
    assert observable_data.nframes is None

    # Check that we can set a observable array, the number of
    # frames gets updated, and the two getter methods are equivalent
    num_frames = 13
    observable_data.kinetic_energy = np.random.random(num_frames)
    assert observable_data.nframes == num_frames
    assert np.array_equal(
        observable_data["kinetic_energy"], observable_data.kinetic_energy
    )

    # Check that we can set another observable array, the number of
    # frames stays correct, and the two getter methods are equivalent
    observable_data["potential_energy"] = np.random.random(num_frames)
    assert observable_data.nframes == num_frames
    assert np.array_equal(
        observable_data["potential_energy"], observable_data.potential_energy
    )

    # Check that we can set an observable array back to None, the number of
    # frames stays correct, and the two getter methods are equivalent
    observable_data.potential_energy = None
    assert observable_data.nframes == num_frames
    assert observable_data.potential_energy is None
    assert np.array_equal(
        observable_data["potential_energy"], observable_data.potential_energy
    )

    # Check that it only accepts one-dimensional arrays
    with pytest.raises(pv_error.InputError):
        # 0 dimensions
        observable_data.potential_energy = np.array(0)
    with pytest.raises(pv_error.InputError):
        # 2 dimensions (reducible)
        observable_data.potential_energy = np.random.random((1, 2))
    with pytest.raises(pv_error.InputError):
        # 2 dimensions (non-reducible)
        observable_data.potential_energy = np.random.random((3, 2))
    with pytest.raises(pv_error.InputError):
        # 3 dimensions
        observable_data.potential_energy = np.random.random((5, 2, 2))

    # Check that random set / get strings raise an error
    key = "anything"
    assert key not in ObservableData.observables()
    with pytest.raises(KeyError):
        observable_data[key] = np.random.random(num_frames)
    with pytest.raises(KeyError):
        _ = observable_data[key]

    # Check that the number of frames wasn't changed
    assert observable_data.nframes == num_frames

    # Check that different number of frames yields warning
    # and `None` number of frames
    with pytest.warns(UserWarning):
        observable_data.potential_energy = np.random.random(num_frames + 1)
    assert observable_data.nframes is None

    # Check that number of frames can be fixed
    observable_data.potential_energy = np.random.random(num_frames)
    assert observable_data.nframes == num_frames

    # Check that all observables can be read in two ways
    for observable in ObservableData.observables():
        observable_data[observable] = np.random.random(num_frames)
    assert np.array_equal(
        observable_data.kinetic_energy, observable_data["kinetic_energy"]
    )
    assert np.array_equal(
        observable_data.potential_energy, observable_data["potential_energy"]
    )
    assert np.array_equal(observable_data.total_energy, observable_data["total_energy"])
    assert np.array_equal(observable_data.volume, observable_data["volume"])
    assert np.array_equal(observable_data.pressure, observable_data["pressure"])
    assert np.array_equal(observable_data.temperature, observable_data["temperature"])
    assert np.array_equal(
        observable_data.constant_of_motion, observable_data["constant_of_motion"]
    )
    assert np.array_equal(
        observable_data.number_of_species, observable_data["number_of_species"]
    )

    # Check that for the number of species, 2D arrays are allowed and don't mess up the frame number
    observable_data.number_of_species = np.random.random((num_frames, 3))
    assert observable_data.nframes == num_frames

    # Check that for the number of species, a 1D array is transformed into 2D array
    temporary_number_of_species = np.random.random((num_frames, 1))
    observable_data.number_of_species = temporary_number_of_species.flatten()
    assert np.array_equal(
        observable_data.number_of_species, temporary_number_of_species
    )
