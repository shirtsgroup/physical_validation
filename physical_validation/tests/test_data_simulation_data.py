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
This file contains tests for the `physical_validation.data.simulation_data` module.
"""
import numpy as np
import pytest

from ..data import (
    EnsembleData,
    FlatfileParser,
    ObservableData,
    SimulationData,
    SystemData,
    TrajectoryData,
)
from ..util import error as pv_error
from .test_system_database import database


class TestEquality:
    r"""
    Tests the equality function of simulation data and its underlying objects
    """

    @staticmethod
    def test_empty() -> None:
        r"""
        Tests that empty objects are equal
        """
        for class_name in [SimulationData, SystemData, ObservableData, TrajectoryData]:
            print(f"Testing empty {class_name.__name__} objects")
            object1 = class_name()
            object2 = class_name()
            assert object1 == object2

        object1 = EnsembleData("NVE", natoms=10, volume=1)
        object2 = EnsembleData("NVE", natoms=10, volume=1)
        assert object1 == object2

    @staticmethod
    def test_object_with_observables() -> None:
        r"""
        Tests that SimulationData objects with ObservableData are behaving as expected
        """
        system_name = "Water900"
        simulation_id_1 = "NPT-lowT-lowP"
        simulation_id_2 = "NPT-highT-lowP"
        system = database.system(system_name)
        parser = FlatfileParser()
        simulation_data_1 = parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble(simulation_id_1),
            system=system.system_data,
            kinetic_ene_file=system.observable_flat_file(
                simulation_id_1, "kinetic_energy"
            ),
            potential_ene_file=system.observable_flat_file(
                simulation_id_1, "potential_energy"
            ),
            total_ene_file=system.observable_flat_file(simulation_id_1, "total_energy"),
            volume_file=system.observable_flat_file(simulation_id_1, "volume"),
        )
        simulation_data_1_copy = parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble(simulation_id_1),
            system=system.system_data,
            kinetic_ene_file=system.observable_flat_file(
                simulation_id_1, "kinetic_energy"
            ),
            potential_ene_file=system.observable_flat_file(
                simulation_id_1, "potential_energy"
            ),
            total_ene_file=system.observable_flat_file(simulation_id_1, "total_energy"),
            volume_file=system.observable_flat_file(simulation_id_1, "volume"),
        )
        simulation_data_2 = parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble(simulation_id_2),
            system=system.system_data,
            kinetic_ene_file=system.observable_flat_file(
                simulation_id_2, "kinetic_energy"
            ),
            potential_ene_file=system.observable_flat_file(
                simulation_id_2, "potential_energy"
            ),
            total_ene_file=system.observable_flat_file(simulation_id_2, "total_energy"),
            volume_file=system.observable_flat_file(simulation_id_2, "volume"),
        )

        assert simulation_data_1 == simulation_data_1_copy
        assert simulation_data_1 != simulation_data_2

    @staticmethod
    def test_object_with_trajectory() -> None:
        r"""
        Tests that SimulationData objects with TrajectoryData are behaving as expected
        """
        system_name = "Octanol2"
        system = database.system(system_name)
        parser = FlatfileParser()
        simulation_data = parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble("GasPhase"),
            system=system.system_data,
            position_file=system.trajectory_flat_file("GasPhase", "position"),
            velocity_file=system.trajectory_flat_file("GasPhase", "velocity"),
        )
        simulation_data_copy = parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble("GasPhase"),
            system=system.system_data,
            position_file=system.trajectory_flat_file("GasPhase", "position"),
            velocity_file=system.trajectory_flat_file("GasPhase", "velocity"),
        )

        assert simulation_data == simulation_data_copy

        # Change a position entry
        simulation_data_copy.trajectory.position[0][0] += 1
        # Suppressed LGTM warning: See https://github.com/github/codeql/issues/5777
        assert simulation_data != simulation_data_copy  # lgtm [py/redundant-comparison]


class TestInvalidityChecks:
    @staticmethod
    def test_empty_raises_error() -> None:
        simulation_data = SimulationData()
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_units_are_none(
                test_name="dummy_test", argument_name="data"
            )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=True,
            )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_observable_data_is_invalid(
                required_observables=[],
                test_name="dummy_test",
                argument_name="data",
            )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_trajectory_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
            )

    @staticmethod
    @pytest.mark.filterwarnings("ignore:NVE|NVT|NPT|muVT with undefined")
    def test_ensemble_check() -> None:
        simulation_data = SimulationData()
        # NVE isn't tested, so passes always
        simulation_data.ensemble = EnsembleData(ensemble="NVE")
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        # NVT doesn't pass without temperature
        simulation_data.ensemble = EnsembleData(ensemble="NVT")
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        # Adding temperature fixes it
        simulation_data.ensemble = EnsembleData(ensemble="NVT", temperature=300)
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        # NPT doesn't pass without temperature
        simulation_data.ensemble = EnsembleData(ensemble="NPT")
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=False,
                check_mu=True,
            )
        # Adding temperature fixes it if check_pressure is False
        simulation_data.ensemble = EnsembleData(ensemble="NPT", temperature=300)
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=False,
            check_mu=True,
        )
        # But not if pressure is required to be tested
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        # Pressure alone doesn't help
        simulation_data.ensemble = EnsembleData(ensemble="NPT", pressure=1.0)
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        # Temperature and pressure passes
        simulation_data.ensemble = EnsembleData(
            ensemble="NPT", temperature=300, pressure=1.0
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=False,
            check_mu=True,
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        # muVT doesn't pass without temperature
        simulation_data.ensemble = EnsembleData(ensemble="muVT")
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=False,
            )
        # But passes with temperature if mu is not tested
        simulation_data.ensemble = EnsembleData(ensemble="muVT", temperature=300)
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=False,
        )
        # And fails again if mu is required
        simulation_data.ensemble = EnsembleData(ensemble="muVT", temperature=300)
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        # And fails also if mu is specified, but not the temperature
        simulation_data.ensemble = EnsembleData(ensemble="muVT", mu=0.8)
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_ensemble_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_pressure=True,
                check_mu=True,
            )
        # With both specified, it passes with different types of mu input
        simulation_data.ensemble = EnsembleData(
            ensemble="muVT", mu=0.8, temperature=300
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        simulation_data.ensemble = EnsembleData(
            ensemble="muVT", mu=[0.8], temperature=300
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        simulation_data.ensemble = EnsembleData(
            ensemble="muVT", mu=np.array([0.8]), temperature=300
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        simulation_data.ensemble = EnsembleData(
            ensemble="muVT", mu=[0.8, 0.5], temperature=300
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )
        simulation_data.ensemble = EnsembleData(
            ensemble="muVT", mu=np.array([0.8, 0.5]), temperature=300
        )
        simulation_data.raise_if_ensemble_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_pressure=True,
            check_mu=True,
        )

    @staticmethod
    def test_system_check() -> None:
        simulation_data = SimulationData()
        simulation_data.system = SystemData()
        # empty fails
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=True,
            )
        # adding only number of atoms still fails
        simulation_data.system.natoms = 100
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=True,
            )
        # adding number of constraints still fails
        simulation_data.system.nconstraints = 300
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=True,
            )
        # adding number of translational constraints still fails
        simulation_data.system.ndof_reduction_tra = 3
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=True,
            )
        # adding number of rotational constraints passes
        simulation_data.system.ndof_reduction_rot = 3
        simulation_data.raise_if_system_data_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_full_system_data_only=True,
        )
        # but fails on entire system check
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=False,
            )
        # adding masses still fails
        simulation_data.system.mass = np.ones(100)
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=False,
            )
        # adding molecule index still fails
        simulation_data.system.molecule_idx = np.linspace(
            0, 100 * 3, 100, endpoint=False, dtype=int
        )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_system_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
                check_full_system_data_only=False,
            )
        # adding number of constraints per molecule passes
        simulation_data.system.nconstraints_per_molecule = 3 * np.ones(100)
        simulation_data.raise_if_system_data_is_invalid(
            test_name="dummy_test",
            argument_name="data",
            check_full_system_data_only=False,
        )

    @staticmethod
    @pytest.mark.filterwarnings("ignore:ObservableData. Mismatch in number of frames")
    def test_observable_check() -> None:
        simulation_data = SimulationData()
        num_frames = 10
        simulation_data.observables = ObservableData()
        # empty object passes if no observables are required
        simulation_data.raise_if_observable_data_is_invalid(
            required_observables=[],
            test_name="dummy_test",
            argument_name="data",
        )
        # empty object fails if at least one quantity is required
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_observable_data_is_invalid(
                required_observables=["kinetic_energy"],
                test_name="dummy_test",
                argument_name="data",
            )
        # adding the required quantity makes it pass
        simulation_data.observables["kinetic_energy"] = np.random.random(num_frames)
        simulation_data.raise_if_observable_data_is_invalid(
            required_observables=["kinetic_energy"],
            test_name="dummy_test",
            argument_name="data",
        )
        # check fails again if additional observable is requested
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_observable_data_is_invalid(
                required_observables=["kinetic_energy", "potential_energy"],
                test_name="dummy_test",
                argument_name="data",
            )
        # check fails also if wrong number of frames is added
        simulation_data.observables["potential_energy"] = np.random.random(
            num_frames - 1
        )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_observable_data_is_invalid(
                required_observables=["kinetic_energy", "potential_energy"],
                test_name="dummy_test",
                argument_name="data",
            )
        # but passes with identical number of frames
        simulation_data.observables["potential_energy"] = np.random.random(num_frames)
        simulation_data.raise_if_observable_data_is_invalid(
            required_observables=["kinetic_energy", "potential_energy"],
            test_name="dummy_test",
            argument_name="data",
        )

    @staticmethod
    def test_trajectory_check() -> None:
        simulation_data = SimulationData()
        num_atoms = 10
        num_frames = 2
        simulation_data.trajectory = TrajectoryData()
        # empty object fails
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_trajectory_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
            )
        # adding only one trajectory still fails
        simulation_data.trajectory.position = np.random.random(
            (num_frames, num_atoms, 3)
        )
        with pytest.raises(pv_error.InputError):
            simulation_data.raise_if_trajectory_data_is_invalid(
                test_name="dummy_test",
                argument_name="data",
            )
        # adding correct velocity trajectory passes
        simulation_data.trajectory.velocity = np.random.random(
            (num_frames, num_atoms, 3)
        )
        simulation_data.raise_if_trajectory_data_is_invalid(
            test_name="dummy_test",
            argument_name="data",
        )
