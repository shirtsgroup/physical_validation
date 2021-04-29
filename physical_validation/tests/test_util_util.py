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
This file contains tests for the `physical_validation.util.util` module.
"""
import numpy as np

from ..util.util import array_equal_shape_and_close


class TestArrayEqualShapeAndClose:
    @staticmethod
    def test_none_is_true() -> None:
        assert array_equal_shape_and_close(None, None)

    @staticmethod
    def test_none_with_array_is_false() -> None:
        assert not array_equal_shape_and_close(None, np.array(None))
        assert not array_equal_shape_and_close(None, np.array(0))
        assert not array_equal_shape_and_close(None, np.zeros(0))
        assert not array_equal_shape_and_close(np.array(None), None)
        assert not array_equal_shape_and_close(np.array(0), None)
        assert not array_equal_shape_and_close(np.zeros(0), None)

    @staticmethod
    def test_shape_and_content() -> None:
        assert array_equal_shape_and_close(np.ones((4, 5)), np.ones((4, 5)))
        assert not array_equal_shape_and_close(np.ones((4, 5)), np.ones((5, 4)))
        assert not array_equal_shape_and_close(np.ones(0), np.ones(1))
        assert array_equal_shape_and_close(np.ones(0), np.zeros(0))
