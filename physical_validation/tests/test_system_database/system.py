###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
The System object describes a system present in the database. It allows
users access to all relevant data used by the physical validation tests
(units, information about the system, information about the ensemble,
and energy and coordinate trajectories).
"""
from typing import Dict, KeysView, Optional, List

import numpy as np

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
        observable_as_array: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
        trajectory_flat_file: Optional[Dict[str, Dict[str, str]]] = None,
        trajectory_as_array: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
    ) -> None:
        if observable_flat_file is None:
            observable_flat_file = {}
        if observable_as_array is None:
            observable_as_array = {}
        if trajectory_flat_file is None:
            trajectory_flat_file = {}
        if trajectory_as_array is None:
            trajectory_as_array = {}

        self.__units = units
        self.__system_data = system_data
        self.__ensemble = ensemble
        self.__description = description
        if simulation_keys == "ensemble":
            self.__simulations = ensemble.keys()
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
        if observable_as_array:
            assert observable_as_array.keys() == self.__simulations

        self.__observable_flat_file = observable_flat_file
        self.__observable_as_array = observable_as_array
        self.__trajectory_flat_file = trajectory_flat_file
        self.__trajectory_as_array = trajectory_as_array

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

    def observable_as_array(
        self, simulation_key: str, quantity: str
    ) -> Optional[np.ndarray]:
        if (
            self.__observable_as_array is None
            or quantity not in self.__observable_as_array[simulation_key]
        ):
            return None
        return self.__observable_as_array[simulation_key][quantity]

    def trajectory_flat_file(self, simulation_key: str, quantity: str) -> Optional[str]:
        if (
            self.__trajectory_flat_file is None
            or simulation_key not in self.__trajectory_flat_file
            or quantity not in self.__trajectory_flat_file[simulation_key]
        ):
            return None
        return self.__trajectory_flat_file[simulation_key][quantity]

    def trajectory_as_array(
        self, simulation_key: str, quantity: str
    ) -> Optional[np.ndarray]:
        if (
            self.__trajectory_as_array is None
            or simulation_key not in self.__trajectory_as_array
            or quantity not in self.__trajectory_as_array[simulation_key]
        ):
            return None
        return self.__trajectory_as_array[simulation_key][quantity]
