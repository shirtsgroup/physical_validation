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
This file contains tests for the `physical_validation.data.lammps_parser` module.
"""
import physical_validation as pv

from .test_data_flatfile_parser import TestFlatFileParser
from .test_system_database import database


class TestLammpsParser:
    r"""
    Test fixture for the LAMMPS parser
    """

    @staticmethod
    def test_read_lammps() -> None:
        system_name = "Water300"
        flat_file_parser = pv.data.FlatfileParser()
        system = database.system(system_name)
        lammps_files = database.lammps_files(system_name)
        parser = pv.data.LammpsParser()
        ensemble = pv.data.EnsembleData(
            ensemble="NVT", natoms=100 * 3, volume=20 ** 3, temperature=300
        )
        lammps_simulation_data = parser.get_simulation_data(
            ensemble=ensemble,
            in_file=lammps_files["in"],
            log_file=lammps_files["log"],
            data_file=lammps_files["data"],
            dump_file=lammps_files["dump"],
        )

        simulation_data_flat_full = TestFlatFileParser.get_flat_file_simulation_data(
            parser=flat_file_parser,
            system=system,
            simulation_id="NVT",
        )

        assert lammps_simulation_data == simulation_data_flat_full
