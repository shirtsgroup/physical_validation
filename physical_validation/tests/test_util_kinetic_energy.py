###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
This file contains tests for the `physical_validation.util.kinetic_energy` module.
"""
import itertools

import pytest

from ..util import error as pv_error
from ..util.kinetic_energy import is_close, temperature


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
        assert is_close(1 + 9.99e-10, 1, rel_tol=1e-9)
        assert not is_close(1 + 1.01e-9, 1, rel_tol=1e-9)
        assert is_close(1e3 + 9.99e-10, 1e3, rel_tol=1e-12)
        assert not is_close(1e3 + 1.01e-9, 1e3, rel_tol=1e-12)

    @staticmethod
    def test_absolute_tolerance_works():
        assert is_close(1e-8, 1.099e-8, abs_tol=1e-9)
        assert not is_close(1e-8, 1.101e-8, abs_tol=1e-9)

    @staticmethod
    def test_small_numbers_work():
        assert is_close(1e-15, 0)
        assert is_close(-9.9e-10, 0)
        assert not is_close(2e-9, 0)
        assert not is_close(-1.1e-9, 0)


class TestTemperature:
    r"""
    Tests the `temperature(kin, ndof, kb)` function
    """

    @staticmethod
    def test_zero_works():
        assert temperature(kin=0, ndof=457.4) == 0
        assert temperature(kin=0, ndof=457.4, kb=8.314e-3) == 0

    @staticmethod
    def test_zero_values_throw():
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=0)
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=0, kb=8.314e-3)
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=6234.9, kb=0)

    @staticmethod
    def test_negative_values_throw():
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=-1e-12)
            temperature(kin=234.1, ndof=-1e-5)
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=-13.6, kb=8.314e-3)
        with pytest.raises(pv_error.InputError):
            temperature(kin=234.1, ndof=6234.9, kb=-8.314e-3)

    @staticmethod
    def test_default_kb():
        # Changing the default would likely yield some ugly surprises
        # TODO: Remove the default - if we support different units,
        #       we should force users to be specific!
        assert temperature(kin=12.634, ndof=76.8, kb=8.314e-3) == temperature(
            kin=12.634, ndof=76.8
        )

    @staticmethod
    def test_formula():
        kinetic_energy = [3468.7, 1.60923e-8, 5.7923e15]
        num_dof = [1, 0.002, 9.6234e13]
        kb_values = [8.314e-3, 1, 1.632e2]
        for [k, n, kb] in itertools.product(kinetic_energy, num_dof, kb_values):
            assert temperature(k, n, kb) == 2 * k / (n * kb)
