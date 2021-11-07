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
The System object describes a system present in the database. It allows
users access to all relevant data used by the physical validation tests
(units, information about the system, information about the ensemble,
and energy and coordinate trajectories).
"""
from typing import Dict, KeysView, List, Optional

import physical_validation as pv


class System:
    def __init__(
        self,
        units: pv.data.UnitData,
        system_data: pv.data.SystemData,
        ensemble: Dict[str, pv.data.EnsembleData],
        description: str,
        simulation_keys: str,
        time_step: Optional[List[float]] = None,
        observable_flat_file: Optional[Dict[str, Dict[str, str]]] = None,
        trajectory_flat_file: Optional[Dict[str, Dict[str, str]]] = None,
    ) -> None:
        if observable_flat_file is None:
            observable_flat_file = {}
        if trajectory_flat_file is None:
            trajectory_flat_file = {}

        self.__units = units
        self.__system_data = system_data
        self.__ensemble = ensemble
        self.__description = description
        if simulation_keys == "ensemble":
            self.__simulations = ensemble.keys()
            if time_step is not None:
                if len(time_step) > 1:
                    raise ValueError(
                        'simulation_keys == "ensemble" can have at most 1 `time_step` input.'
                    )
                if len(time_step) == 1:
                    self.__time_step = {key: time_step[0] for key in self.__simulations}
        elif simulation_keys == "time step":
            if time_step is None:
                raise ValueError(
                    'simulation_keys == "time step" requires `time_step` input.'
                )
            self.__time_step = {str(dt): dt for dt in time_step}
            self.__simulations = self.__time_step.keys()
        else:
            raise NotImplementedError("Unknown simulation keys.")

        if observable_flat_file:
            assert observable_flat_file.keys() == self.__simulations

        self.__observable_flat_file = observable_flat_file
        self.__trajectory_flat_file = trajectory_flat_file

    @property
    def units(self) -> pv.data.UnitData:
        return self.__units

    @property
    def system_data(self) -> pv.data.SystemData:
        return self.__system_data

    @property
    def simulations(self) -> KeysView[str]:
        return self.__simulations

    @property
    def description(self) -> str:
        return self.__description

    def ensemble(self, simulation_key: str) -> pv.data.EnsembleData:
        return self.__ensemble[simulation_key]

    def time_step(self, simulation_key: str) -> float:
        return self.__time_step[simulation_key]

    def observable_flat_file(self, simulation_key: str, quantity: str) -> Optional[str]:
        if (
            self.__observable_flat_file is None
            or quantity not in self.__observable_flat_file[simulation_key]
        ):
            return None
        return self.__observable_flat_file[simulation_key][quantity]

    def trajectory_flat_file(self, simulation_key: str, quantity: str) -> Optional[str]:
        if (
            self.__trajectory_flat_file is None
            or simulation_key not in self.__trajectory_flat_file
            or quantity not in self.__trajectory_flat_file[simulation_key]
        ):
            return None
        return self.__trajectory_flat_file[simulation_key][quantity]
