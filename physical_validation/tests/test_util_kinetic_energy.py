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
This file contains tests for the `physical_validation.util.kinetic_energy` module.
"""
from ..util.kinetic_energy import is_close


class TestIsClose:
    r"""
    Tests the `is_close(a, b, rel_tol, abs_tol)` function
    """

    @staticmethod
    def test_identical_values_are_close():
        assert is_close(0, 0)
        assert is_close(-15, -15)

    @staticmethod
    def test_significantly_different_values_are_not_close():
        assert not is_close(1, 0)
        assert not is_close(1, 0)
        assert not is_close(-3, 3)

    @staticmethod
    def test_relative_tolerance_works():
        assert is_close(1 + 9.99e-10, 1, relative_tolerance=1e-9)
        assert not is_close(1 + 1.01e-9, 1, relative_tolerance=1e-9)
        assert is_close(1e3 + 9.99e-10, 1e3, relative_tolerance=1e-12)
        assert not is_close(1e3 + 1.01e-9, 1e3, relative_tolerance=1e-12)

    @staticmethod
    def test_absolute_tolerance_works():
        assert is_close(1e-8, 1.099e-8, absolute_tolerance=1e-9)
        assert not is_close(1e-8, 1.101e-8, absolute_tolerance=1e-9)

    @staticmethod
    def test_small_numbers_work():
        assert is_close(1e-15, 0)
        assert is_close(-9.9e-10, 0)
        assert not is_close(2e-9, 0)
        assert not is_close(-1.1e-9, 0)
