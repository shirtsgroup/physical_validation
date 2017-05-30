from __future__ import print_function
from __future__ import division

import numpy as np


def is_num_ndarray(array, n=None, min_n=None, max_n=None):
    assert isinstance(array, np.ndarray)
    assert np.issubdtype(array.dtype, np.number)
    if n is not None:
        assert array.ndim == n
    if min_n is not None:
        assert array.ndim >= min_n
    if max_n is not None:
        assert array.ndim <= max_n
    return True


def have_identical_shape(array1, array2):
    assert isinstance(array1, np.ndarray)
    assert isinstance(array2, np.ndarray)

    assert array1.shape == array2.shape

    return True


def has_shape(array, shape):
    assert isinstance(array, np.ndarray)

    assert array.shape == shape

    return True


def equipartition(positions, velocities, masses, molec_idx=None, molec_nbonds=None):

    # masses: 1 x natoms
    masses = np.array(masses)
    is_num_ndarray(masses, 1)
    natoms = masses.shape[0]

    # positions / velocities: nframes x natoms x 3
    positions = np.array(positions)
    velocities = np.array(velocities)
    is_num_ndarray(positions, min_n=2, max_n=3)
    is_num_ndarray(velocities, min_n=2, max_n=3)
    have_identical_shape(positions, velocities)

    if positions.ndim == 2:
        nframes = 1
        has_shape(positions, (natoms, 3))
    else:
        nframes = positions.shape[0]
        has_shape(positions, (nframes, natoms, 3))

    if molec_idx is not None:
        molec_idx = np.array(molec_idx)
        is_num_ndarray(molec_idx, 1)
        nmolecs = molec_idx.shape[0]
        assert nmolecs <= natoms
        assert molec_idx[-1] < natoms

    if molec_nbonds is not None:
        molec_nbonds = np.array(molec_nbonds)
        is_num_ndarray(molec_nbonds, 1)
        assert molec_nbonds.shape[0] == nmolecs

    return (positions, velocities, masses,
            molec_idx, molec_nbonds,
            natoms, nmolecs, nframes)
