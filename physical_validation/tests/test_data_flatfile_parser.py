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
This file contains tests for the `physical_validation.data.gromacs_parser` module.
"""

import numpy as np

from ..data import FlatfileParser, ObservableData, SimulationData, TrajectoryData
from .test_system_database import database


class TestFlatFileParser:
    r"""
    Test fixture for the flat file parser
    """

    @staticmethod
    def get_flat_file_simulation_data(
        parser: FlatfileParser, system: database.System, simulation_id: str
    ) -> SimulationData:
        r"""
        Returns a SimulationData object created with the FlatFileParser

        Parameters
        ----------
        parser
            The FlatFileParser object
        system
            A testing System object to read the input from
        simulation_id
            The simulation id to read the input of

        Returns
        -------
        SimulationData object created by the flat file parser
        """
        return parser.get_simulation_data(
            units=system.units,
            ensemble=system.ensemble(simulation_id),
            system=system.system_data,
            dt=system.time_step(simulation_id),
            position_file=system.trajectory_flat_file(simulation_id, "position"),
            velocity_file=system.trajectory_flat_file(simulation_id, "velocity"),
            kinetic_ene_file=system.observable_flat_file(
                simulation_id, "kinetic_energy"
            ),
            potential_ene_file=system.observable_flat_file(
                simulation_id, "potential_energy"
            ),
            total_ene_file=system.observable_flat_file(simulation_id, "total_energy"),
            volume_file=system.observable_flat_file(simulation_id, "volume"),
            temperature_file=system.observable_flat_file(simulation_id, "temperature"),
            pressure_file=system.observable_flat_file(simulation_id, "pressure"),
            const_of_mot_file=system.observable_flat_file(
                simulation_id, "constant_of_motion"
            ),
            number_of_species_file=system.observable_flat_file(
                simulation_id, "number_of_species"
            ),
        )

    @staticmethod
    def get_empty_simulation_data(
        system: database.System, simulation_id: str
    ) -> SimulationData:
        r"""
        Returns a SimulationData object filled with unit, ensemble,
        system and time step information from the testing System object

        Parameters
        ----------
        system
            A testing System object to read the input from
        simulation_id
            The simulation id to read the input of

        Returns
        -------
        SimulationData object without observables or trajectory
        """
        return SimulationData(
            units=system.units,
            ensemble=system.ensemble(simulation_id),
            system=system.system_data,
            dt=system.time_step(simulation_id),
            observables=ObservableData(),
            trajectory=TrajectoryData(),
        )

    @staticmethod
    def test_flat_file_input() -> None:
        r"""
        Tests that the flat file parser returns a SimulationData object that is
        equivalent to creating the object by hand from numpy arrays
        """
        system_name = "Water5"
        simulation_id = "NPT full trajectory"
        flat_file_parser = FlatfileParser()
        system = database.system(system_name)
        simulation_data_flat = TestFlatFileParser.get_flat_file_simulation_data(
            parser=flat_file_parser, system=system, simulation_id=simulation_id
        )
        simulation_data_reference = TestFlatFileParser.get_empty_simulation_data(
            system=system, simulation_id=simulation_id
        )

        # Read observable data from flat file
        for quantity in ObservableData.observables():
            simulation_data_reference.observables[quantity] = np.loadtxt(
                system.observable_flat_file(
                    simulation_key=simulation_id, quantity=quantity
                )
            )
        if simulation_data_reference.observables.number_of_species.ndim == 1:
            simulation_data_reference.observables.number_of_species = (
                simulation_data_reference.observables.number_of_species[:, np.newaxis]
            )

        # Read xyz data from flat file
        for quantity in ["position", "velocity"]:
            trajectory = []
            with open(
                system.trajectory_flat_file(
                    simulation_key=simulation_id, quantity=quantity
                )
            ) as input_file:
                frame = []
                for line in input_file:
                    line = line.strip()
                    if not line:
                        # empty line, frame separator
                        if frame:
                            trajectory.append(frame)
                            frame = []
                        continue
                    line = line.split("#", 1)[0]
                    if not line:
                        # only comment on this line
                        continue
                    line = line.split()
                    assert len(line) == 3
                    frame.append([float(line[0]), float(line[1]), float(line[2])])
                if frame:
                    trajectory.append(frame)
            simulation_data_reference.trajectory[quantity] = np.array(trajectory)

        # ensure flat file parser did the right thing
        assert simulation_data_flat == simulation_data_reference
