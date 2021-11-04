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
from ..data import EnsembleData, ObservableData, SimulationData, UnitData


def ensemble(ensemble_string):
    r"""
    Helper function creating dummy EnsembleData objects for the tests

    Parameters
    ----------
    ensemble_string
        Currently only accepts "NVE"

    Returns
    -------
    An EnsembleData object of ensemble `ensemble_string`
    """
    if ensemble_string == "NVE":
        return EnsembleData(ensemble_string, natoms=10, volume=1.0, energy=5.0)
    raise NotImplementedError(f"Unknown ensemble {ensemble_string}")


class TestUnimplementedEnsemblesThrow:
    r"""
    Tests that unimplemented (and hence untested) ensembles throw.
    If you remove an ensemble from this test, you're expected to add
    some proper tests for it!
    """

    @staticmethod
    def test_ensemble_check_throws():
        for ensemble_string in ["NVE"]:
            simulation_data_1 = SimulationData()
            simulation_data_1.ensemble = ensemble(ensemble_string)
            simulation_data_1.units = UnitData.units("GROMACS")
            simulation_data_2 = SimulationData()
            simulation_data_2.ensemble = ensemble(ensemble_string)
            simulation_data_2.units = UnitData.units("GROMACS")

            with pytest.raises(NotImplementedError):
                pv_ensemble.check(simulation_data_1, simulation_data_2)

    @staticmethod
    def test_interval_estimate_throws():
        for ensemble_string in ["NVE"]:
            simulation_data_1 = SimulationData()
            simulation_data_1.ensemble = ensemble(ensemble_string)
            simulation_data_1.observables = ObservableData()
            simulation_data_1.units = UnitData.units("GROMACS")

            with pytest.raises(NotImplementedError):
                pv_ensemble.estimate_interval(simulation_data_1)
