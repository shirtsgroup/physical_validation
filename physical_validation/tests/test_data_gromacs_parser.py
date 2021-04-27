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

from typing import Dict

import physical_validation as pv

from .test_system_database import database


class TestGromacsParser:
    r"""
    Test fixture for the GROMACS parser
    """

    @staticmethod
    def get_flat_file_simulation_data(
        parser: pv.data.FlatfileParser, system: database.System, simulation_id: str
    ) -> pv.data.SimulationData:
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
        )

    @staticmethod
    def get_gromacs_simulation_data(
        parser: pv.data.GromacsParser,
        gromacs_files: Dict[str, str],
        use_full_trajectory: bool,
    ) -> pv.data.SimulationData:
        return parser.get_simulation_data(
            mdp=gromacs_files["parameters"],
            top=gromacs_files["topology"],
            edr=gromacs_files["energy"],
            gro=gromacs_files["final configuration"]
            if not use_full_trajectory
            else None,
            trr=gromacs_files["trajectory"] if use_full_trajectory else None,
        )

    @staticmethod
    def test_gromacs_and_flat_file_parsers_are_equivalent() -> None:
        r"""
        Test comparing simulation data objects created with the GROMACS
        and the flat file parser, respectively
        """
        system_name = "Water5"
        flat_file_parser = pv.data.FlatfileParser()
        system = database.system(system_name)
        gromacs_parser = pv.data.GromacsParser(
            "/opt/anaconda3/envs/pv-test/bin.AVX2_256/gmx"
        )
        gromacs_files = database.gromacs_files(system_name)

        simulation_data_flat_full = TestGromacsParser.get_flat_file_simulation_data(
            parser=flat_file_parser, system=system, simulation_id="full trajectory"
        )
        simulation_data_gromacs_full = TestGromacsParser.get_gromacs_simulation_data(
            parser=gromacs_parser, gromacs_files=gromacs_files, use_full_trajectory=True
        )
        assert simulation_data_gromacs_full == simulation_data_flat_full

        simulation_data_flat_last_frame = (
            TestGromacsParser.get_flat_file_simulation_data(
                parser=flat_file_parser, system=system, simulation_id="last frame only"
            )
        )
        simulation_data_gromacs_last_frame = (
            TestGromacsParser.get_gromacs_simulation_data(
                parser=gromacs_parser,
                gromacs_files=gromacs_files,
                use_full_trajectory=False,
            )
        )
        assert simulation_data_gromacs_last_frame == simulation_data_flat_last_frame
