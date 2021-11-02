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
This file contains tests for the `physical_validation.util.ensemble` module.
"""
import numpy as np

from ..util.ensemble import chemical_potential_energy


class TestChemicalPotentialEnergy:
    r"""
    Tests the `chemical_potential_energy` function
    """

    @staticmethod
    def test_1d():
        mu = np.array([0.4])
        num_species = np.array([3, 6, 7])
        expected_result = mu * num_species
        assert np.all(chemical_potential_energy(mu, num_species) == expected_result)

    @staticmethod
    def test_2d():
        mu = np.array([0.3, 0.8, 0.7])
        num_species = np.array([[2, 4, 3], [9, 1, 3], [5, 6, 3], [3, 5, 3], [4, 2, 5]])
        expected_result = np.sum(mu * num_species, axis=1)
        assert np.all(chemical_potential_energy(mu, num_species) == expected_result)
