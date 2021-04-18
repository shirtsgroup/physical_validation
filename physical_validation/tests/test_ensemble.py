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
This file contains tests for the `physical_validation.ensemble` module.
"""
import pytest

from .. import ensemble as pv_ensemble
from ..data import EnsembleData, ObservableData, SimulationData
from ..util import error as pv_error


class TestUnimplementedEnsemblesThrow:
    r"""
    Tests that unimplemented (and hence untested) ensembles throw.
    If you remove an ensemble from this test, you're expected to add
    some proper tests for it!
    """

    @staticmethod
    def test_ensemble_check_throws():
        for ensemble in ["NVE", "muVT"]:
            simulation_data_1 = SimulationData()
            simulation_data_1.ensemble = EnsembleData(ensemble)
            simulation_data_2 = SimulationData()
            simulation_data_2.ensemble = EnsembleData(ensemble)

            with pytest.raises(pv_error.InputError):
                pv_ensemble.check(simulation_data_1, simulation_data_2)

    @staticmethod
    def test_interval_estimate_throws():
        for ensemble in ["NVE", "muVT"]:
            simulation_data_1 = SimulationData()
            simulation_data_1.ensemble = EnsembleData(ensemble)
            simulation_data_1.observables = ObservableData()

            with pytest.raises(NotImplementedError):
                pv_ensemble.estimate_interval(simulation_data_1)
