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
This file contains tests for the `physical_validation.integrator` module.
"""
import pytest

from .. import integrator
from ..data import ObservableData, SimulationData
from ..util import error as pv_error


class TestInvalidInputThrows:
    r"""
    Tests that invalid input throws.
    """

    @staticmethod
    def test_non_simulation_data_throws():
        r"""
        Tests that non-SimulationData objects in the list get detected.
        """
        observable_data = ObservableData(constant_of_motion=[1])
        simulation_data_list = [
            SimulationData(dt=0.002, observables=observable_data),
            SimulationData(dt=0.001, observables=observable_data),
            [],
        ]
        with pytest.raises(pv_error.InputError):
            integrator.convergence(simulation_data_list)
        simulation_data_list[-1] = None
        with pytest.raises(pv_error.InputError):
            integrator.convergence(simulation_data_list)
        simulation_data_list[-1] = 0
        with pytest.raises(pv_error.InputError):
            integrator.convergence(simulation_data_list)

    @staticmethod
    def test_invalid_time_step_throws():
        r"""
        Tests that invalid time steps in the SimulationData object get detected.
        """
        simulation_data = SimulationData()
        simulation_data.dt = 0
        with pytest.raises(pv_error.InputError):
            integrator.convergence([simulation_data])
        simulation_data.dt = -0.1
        with pytest.raises(pv_error.InputError):
            integrator.convergence([simulation_data])
        simulation_data.dt = -1
        with pytest.raises(pv_error.InputError):
            integrator.convergence([simulation_data])

    @staticmethod
    def test_identical_time_step_throws():
        r"""
        Tests that identical time steps in the SimulationData object list get detected.
        """
        observable_data = ObservableData(constant_of_motion=[1])
        simulation_data_list = [
            SimulationData(dt=0.1, observables=observable_data),
            SimulationData(dt=0.2, observables=observable_data),
            SimulationData(dt=0.1, observables=observable_data),
        ]
        with pytest.raises(pv_error.InputError):
            integrator.convergence(simulation_data_list)
