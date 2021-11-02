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
flatfile_parser.py
"""
from typing import List, Optional

from . import (
    EnsembleData,
    ObservableData,
    SimulationData,
    SystemData,
    TrajectoryData,
    UnitData,
    parser,
)


class FlatfileParser(parser.Parser):
    """
    FlatfileParser
    """

    def __init__(self):
        super(FlatfileParser, self).__init__()

    def get_simulation_data(
        self,
        units: Optional[UnitData] = None,
        ensemble: Optional[EnsembleData] = None,
        system: Optional[SystemData] = None,
        dt: Optional[float] = None,
        position_file: Optional[str] = None,
        velocity_file: Optional[str] = None,
        kinetic_ene_file: Optional[str] = None,
        potential_ene_file: Optional[str] = None,
        total_ene_file: Optional[str] = None,
        volume_file: Optional[str] = None,
        pressure_file: Optional[str] = None,
        temperature_file: Optional[str] = None,
        const_of_mot_file: Optional[str] = None,
        number_of_species_file: Optional[str] = None,
    ) -> SimulationData:
        r"""Read simulation data from flat files

        Returns a SimulationData object created from (optionally) provided UnitData, EnsembleData
        and SystemData, as well as TrajectoryData and ObservableData objects created from flat
        files. The files are expected to be in one of the following formats:

        * xyz-format
          trajectory files (position_file, velocity_file)
          - three numbers per line, separated by white space
          - frames delimited by a completely blank line
          - any character after (and including) a '#' are ignored
        * 1d-format
          all other files
          - one number per line
          - any character after (and including) a '#' are ignored

        Parameters
        ----------
        units: UnitData, optional
            A UnitData object representing the units used in the simulation
        ensemble: EnsembleData, optional
            A EnsembleData object representing the ensemble the simulation has been performed in
        system: SystemData, optional
            A SystemData object representing the atoms and molecules in the system
        dt: float, optional
            The time step used in the simulation
        position_file: str, optional
            Path to a file in xyz-format containing the position trajectory
        velocity_file: str, optional
            Path to a file in xyz-format containing the velocity trajectory
        kinetic_ene_file: str, optional
            Path to a file in 1d-format containing the kinetic energy trajectory
        potential_ene_file: str, optional
            Path to a file in 1d-format containing the potential energy trajectory
        total_ene_file: str, optional
            Path to a file in 1d-format containing the total energy trajectory
        volume_file: str, optional
            Path to a file in 1d-format containing the volume trajectory
        pressure_file: str, optional
            Path to a file in 1d-format containing the pressure trajectory
        temperature_file: str, optional
            Path to a file in 1d-format containing the temperature trajectory
        const_of_mot_file: str, optional
            Path to a file in 1d-format containing the constant of motion trajectory
        number_of_species_file: str, optional
            Path to a file in nd-format containing the number of species trajectory

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the provided ensemble and
            system objects as well as the trajectory data found in the
            edr and trr / gro files.

        """

        trj_dict = {"position": position_file, "velocity": velocity_file}

        if any(trj_dict.values()):
            trajectory = TrajectoryData()
            for key, filename in trj_dict.items():
                if filename is None:
                    continue
                trajectory[key] = self.__read_xyz(filename)
        else:
            trajectory = None

        obs_dict = {
            "kinetic_energy": kinetic_ene_file,
            "potential_energy": potential_ene_file,
            "total_energy": total_ene_file,
            "volume": volume_file,
            "pressure": pressure_file,
            "temperature": temperature_file,
            "constant_of_motion": const_of_mot_file,
            "number_of_species": number_of_species_file,
        }

        if any(obs_dict.values()):
            observables = ObservableData()
            for key, filename in obs_dict.items():
                if filename is None:
                    continue
                observables[key] = (
                    self.__read_1d(filename)
                    if key != "number_of_species"
                    else self.__read_nd(filename)
                )
        else:
            observables = None

        result = SimulationData(
            units=units,
            dt=dt,
            system=system,
            ensemble=ensemble,
            observables=observables,
            trajectory=trajectory,
        )

        return result

    @staticmethod
    def __read_xyz(filename: str) -> List[List[List[float]]]:
        result = []
        with open(filename) as f:
            frame = []
            for line in f:
                line = line.strip()
                if not line:
                    # blank line
                    if frame:
                        result.append(frame)
                        frame = []
                    continue
                line = line.split("#", maxsplit=1)[0].strip()
                if not line:
                    # only comment on this line
                    continue
                xyz = line.split()[0:3]
                frame.append([float(n) for n in xyz])
            if frame:
                result.append(frame)
        return result

    @staticmethod
    def __read_1d(filename: str) -> List[float]:
        result = []
        with open(filename) as f:
            for line in f:
                line = line.split("#", maxsplit=1)[0].strip()
                if not line:
                    # blank or comment-only line
                    continue
                result.append(float(line.strip()))
        return result

    @staticmethod
    def __read_nd(filename: str) -> List[List[float]]:
        result = []
        with open(filename) as f:
            for line in f:
                line = line.split("#", maxsplit=1)[0].strip()
                if not line:
                    # blank or comment-only line
                    continue
                result.append([float(entry.strip()) for entry in line.split()])
        return result
