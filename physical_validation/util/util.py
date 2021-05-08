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
Miscellaneous utility functions
"""
from typing import Optional

import numpy as np


def array_equal_shape_and_close(
    array1: Optional[np.ndarray], array2: Optional[np.ndarray]
) -> bool:
    r"""
    Tests whether two arrays have the same shape and all values are close

    Parameters
    ----------
    array1
    array2

    Returns
    -------
    True   if the two arrays have the same shape and all values are close,
    False  otherwise
    """
    if array1 is None and array2 is None:
        return True
    if array1 is None or array2 is None:
        return False
    if array1.shape != array2.shape:
        return False
    return np.allclose(array1, array2, rtol=1e-12, atol=1e-12)
