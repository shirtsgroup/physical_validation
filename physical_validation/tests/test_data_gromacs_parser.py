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
import os
from typing import Dict

import pytest

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
    @pytest.mark.parametrize("ensemble", ["NVT", "NPT"])
    def test_gromacs_and_flat_file_parsers_are_equivalent(ensemble: str) -> None:
        r"""
        Test comparing simulation data objects created with the GROMACS
        and the flat file parser, respectively
        """
        system_name = "Water5"
        flat_file_parser = pv.data.FlatfileParser()
        system = database.system(system_name)
        gromacs_parser = pv.data.GromacsParser()
        gromacs_files = database.gromacs_files(system_name)[ensemble]

        simulation_data_flat_full = TestFlatFileParser.get_flat_file_simulation_data(
            parser=flat_file_parser,
            system=system,
            simulation_id=f"{ensemble} full trajectory",
        )
        # GROMACS doesn't write out the number of species
        simulation_data_flat_full.observables.number_of_species = None
        simulation_data_gromacs_full = TestGromacsParser.get_gromacs_simulation_data(
            parser=gromacs_parser,
            gromacs_files=gromacs_files,
            use_full_trajectory=True,
        )
        assert simulation_data_gromacs_full == simulation_data_flat_full

        simulation_data_flat_last_frame = (
            TestFlatFileParser.get_flat_file_simulation_data(
                parser=flat_file_parser,
                system=system,
                simulation_id=f"{ensemble} last frame only",
            )
        )
        # GROMACS doesn't write out the number of species
        simulation_data_flat_last_frame.observables.number_of_species = None
        simulation_data_gromacs_last_frame = (
            TestGromacsParser.get_gromacs_simulation_data(
                parser=gromacs_parser,
                gromacs_files=gromacs_files,
                use_full_trajectory=False,
            )
        )
        assert simulation_data_gromacs_last_frame == simulation_data_flat_last_frame

    @staticmethod
    def test_gromacs_topology_exception() -> None:
        r"""
        Check that the GROMACS parser raises an exception when a topology
        include cannot be found.
        """
        system_name = "Octanol512"
        parser = pv.data.GromacsParser()
        gromacs_files = database.gromacs_files(system_name)["GasPhase"]
        with pytest.raises(IOError):
            parser.get_simulation_data(
                mdp=gromacs_files["parameters"], top=gromacs_files["topology"]
            )

    @staticmethod
    @pytest.mark.filterwarnings("ignore:NVT with undefined volume")
    def test_gromacs_topology_with_bonds() -> None:
        r"""
        Check that GROMACS parser reads a system with bonds and angle
        terms correctly.
        """
        system_name = "Octanol512"
        system = database.system(system_name)
        gromacs_files = database.gromacs_files(system_name)["GasPhase"]

        parser = pv.data.GromacsParser(includepath=gromacs_files["include_path"])
        simulation_data = parser.get_simulation_data(
            mdp=gromacs_files["parameters"], top=gromacs_files["topology"]
        )
        assert simulation_data.system == system.system_data

    @staticmethod
    def run_gromacs_simulation(
        mdp_file: str, top_file: str, conf_file: str, deffnm: str
    ) -> None:
        gromacs_interface = pv.util.gromacs_interface.GromacsInterface()
        tpr_file = f"{deffnm}.tpr"
        return_code = gromacs_interface.grompp(
            mdp=mdp_file, top=top_file, gro=conf_file, tpr=tpr_file
        )
        assert return_code == 0

        return_code = gromacs_interface.mdrun(tpr=tpr_file, deffnm=deffnm)
        assert return_code == 0

    @staticmethod
    def clean_up_gromacs_run(deffnm: str) -> None:
        def delete(file: str) -> None:
            try:
                os.remove(file)
            except OSError:
                pass

        delete("mdout.mdp")
        delete(f"{deffnm}.gro")
        delete(f"{deffnm}.edr")
        delete(f"{deffnm}.trr")
        delete(f"{deffnm}.log")
        delete(f"{deffnm}.tpr")
        delete(f"{deffnm}.cpt")

    @staticmethod
    def test_run_gromacs() -> None:
        system_name = "Water5"
        gromacs_files = database.gromacs_files(system_name)["NPT"]
        deffnm = "system"
        TestGromacsParser.run_gromacs_simulation(
            mdp_file=gromacs_files["parameters"],
            top_file=gromacs_files["topology"],
            conf_file=gromacs_files["configuration"],
            deffnm=deffnm,
        )
        TestGromacsParser.clean_up_gromacs_run(deffnm=deffnm)
