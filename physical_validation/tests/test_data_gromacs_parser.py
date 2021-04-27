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

from .test_data_flatfile_parser import TestFlatFileParser
from .test_system_database import database


class TestGromacsParser:
    r"""
    Test fixture for the GROMACS parser
    """

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
        gromacs_parser = pv.data.GromacsParser()
        gromacs_files = database.gromacs_files(system_name)

        simulation_data_flat_full = TestFlatFileParser.get_flat_file_simulation_data(
            parser=flat_file_parser, system=system, simulation_id="full trajectory"
        )
        simulation_data_gromacs_full = TestGromacsParser.get_gromacs_simulation_data(
            parser=gromacs_parser, gromacs_files=gromacs_files, use_full_trajectory=True
        )
        assert simulation_data_gromacs_full == simulation_data_flat_full

        simulation_data_flat_last_frame = (
            TestFlatFileParser.get_flat_file_simulation_data(
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
