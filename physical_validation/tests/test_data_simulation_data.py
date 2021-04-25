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
from ..data import (
    EnsembleData,
    FlatfileParser,
    ObservableData,
    SimulationData,
    SystemData,
    TrajectoryData,
)
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
        assert simulation_data != simulation_data_copy
