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
Data structures carrying simulation data.
"""
import warnings
from typing import Any, List, Optional

import numpy as np

from ..util import error as pv_error
from ..util.util import array_equal_shape_and_close


class ObservableData(object):
    r"""ObservableData: The trajectory of (macroscopic) observables during the simulation

    Stores a number of different observables:

    * kinetic_energy: the kinetic energy of the system,
    * potential_energy: the potential energy of the system,
    * total_energy: the total energy of the system,
    * volume: the volume of the system box,
    * pressure: the pressure of the system,
    * temperature: the temperature of the system,
    * constant_of_motion: a constant of motion of the trajectory.

    The observable trajectories can be accessed either using the getters of an object, as in
        observables.kinetic_energy
    or using the key notation, as in
        observables['kinetic_energy']
    """

    @staticmethod
    def observables() -> List[str]:
        return [
            "kinetic_energy",
            "potential_energy",
            "total_energy",
            "volume",
            "pressure",
            "temperature",
            "constant_of_motion",
            "number_of_species",
        ]

    def __init__(
        self,
        kinetic_energy: Any = None,
        potential_energy: Any = None,
        total_energy: Any = None,
        volume: Any = None,
        pressure: Any = None,
        temperature: Any = None,
        constant_of_motion: Any = None,
        number_of_species: Any = None,
    ):
        self.__kinetic_energy = None
        self.__potential_energy = None
        self.__total_energy = None
        self.__volume = None
        self.__pressure = None
        self.__temperature = None
        self.__constant_of_motion = None
        self.__kinetic_energy_per_molec = None
        self.__number_of_species = None

        self.__getters = {
            "kinetic_energy": ObservableData.kinetic_energy.__get__,
            "potential_energy": ObservableData.potential_energy.__get__,
            "total_energy": ObservableData.total_energy.__get__,
            "volume": ObservableData.volume.__get__,
            "pressure": ObservableData.pressure.__get__,
            "temperature": ObservableData.temperature.__get__,
            "constant_of_motion": ObservableData.constant_of_motion.__get__,
            "number_of_species": ObservableData.number_of_species.__get__,
        }

        self.__setters = {
            "kinetic_energy": ObservableData.kinetic_energy.__set__,
            "potential_energy": ObservableData.potential_energy.__set__,
            "total_energy": ObservableData.total_energy.__set__,
            "volume": ObservableData.volume.__set__,
            "pressure": ObservableData.pressure.__set__,
            "temperature": ObservableData.temperature.__set__,
            "constant_of_motion": ObservableData.constant_of_motion.__set__,
            "number_of_species": ObservableData.number_of_species.__set__,
        }

        # Consistency check
        assert set(self.__getters.keys()) == set(self.__setters.keys())
        assert set(self.__getters.keys()) == set(ObservableData.observables())

        self.kinetic_energy = kinetic_energy
        self.potential_energy = potential_energy
        self.total_energy = total_energy
        self.volume = volume
        self.pressure = pressure
        self.temperature = temperature
        self.constant_of_motion = constant_of_motion
        self.number_of_species = number_of_species

    def __getitem__(self, key: str) -> Optional[np.ndarray]:
        if key not in self.observables():
            raise KeyError
        return self.__getters[key](self)

    def __setitem__(self, key: str, value: Any) -> None:
        if key not in self.observables():
            raise KeyError
        self.__setters[key](self, value)

    def __check_value(self, value: Any, key: str) -> Optional[np.ndarray]:
        if value is None:
            return None
        value = np.array(value)
        if value.ndim != 1 and key != "number_of_species":
            raise pv_error.InputError(key, "Expected 1-dimensional array.")
        if key == "number_of_species" and value.ndim == 1:
            # number of species is a 2D array to allow for multiple species
            value = value[:, np.newaxis]
        if self.nframes is not None and self.nframes != value.shape[0]:
            warnings.warn(
                "ObservableData: Mismatch in number of frames. Setting `nframes = None`."
            )
        return value

    @property
    def kinetic_energy(self) -> Optional[np.ndarray]:
        """Get kinetic_energy"""
        return self.__kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, kinetic_energy: Any) -> None:
        """Set kinetic_energy"""
        self.__kinetic_energy = self.__check_value(kinetic_energy, "kinetic_energy")

    @property
    def potential_energy(self) -> Optional[np.ndarray]:
        """Get potential_energy"""
        return self.__potential_energy

    @potential_energy.setter
    def potential_energy(self, potential_energy: Any) -> None:
        """Set potential_energy"""
        self.__potential_energy = self.__check_value(
            potential_energy, "potential_energy"
        )

    @property
    def total_energy(self) -> Optional[np.ndarray]:
        """Get total_energy"""
        return self.__total_energy

    @total_energy.setter
    def total_energy(self, total_energy: Any) -> None:
        """Set total_energy"""
        self.__total_energy = self.__check_value(total_energy, "total_energy")

    @property
    def volume(self) -> Optional[np.ndarray]:
        """Get volume"""
        return self.__volume

    @volume.setter
    def volume(self, volume: Any) -> None:
        """Set volume"""
        self.__volume = self.__check_value(volume, "volume")

    @property
    def pressure(self) -> Optional[np.ndarray]:
        """Get pressure"""
        return self.__pressure

    @pressure.setter
    def pressure(self, pressure: Any) -> None:
        """Set pressure"""
        self.__pressure = self.__check_value(pressure, "pressure")

    @property
    def temperature(self) -> Optional[np.ndarray]:
        """Get temperature"""
        return self.__temperature

    @temperature.setter
    def temperature(self, temperature: Any) -> None:
        """Set temperature"""
        self.__temperature = self.__check_value(temperature, "temperature")

    @property
    def constant_of_motion(self) -> Optional[np.ndarray]:
        """Get constant_of_motion"""
        return self.__constant_of_motion

    @constant_of_motion.setter
    def constant_of_motion(self, constant_of_motion: Any) -> None:
        """Set constant_of_motion"""
        self.__constant_of_motion = self.__check_value(
            constant_of_motion, "constant_of_motion"
        )

    @property
    def number_of_species(self) -> Optional[np.ndarray]:
        """Get number_of_species"""
        return self.__number_of_species

    @number_of_species.setter
    def number_of_species(self, number_of_species: Any) -> None:
        """Set number_of_species"""
        self.__number_of_species = self.__check_value(
            number_of_species, "number_of_species"
        )

    @property
    def nframes(self) -> Optional[int]:
        """Get number of frames"""
        frames = None
        for observable in ObservableData.observables():
            if self[observable] is not None:
                if frames is not None:
                    if self[observable].shape[0] == frames:
                        continue
                    else:
                        return None
                else:
                    frames = self[observable].shape[0]

        return frames

    @property
    def kinetic_energy_per_molecule(self) -> Optional[np.ndarray]:
        """Get kinetic_energy per molecule - used internally"""
        return self.__kinetic_energy_per_molec

    @kinetic_energy_per_molecule.setter
    def kinetic_energy_per_molecule(self, kinetic_energy: Optional[np.ndarray]) -> None:
        """Set kinetic_energy per molecule - used internally"""
        self.__kinetic_energy_per_molec = kinetic_energy

    def __eq__(self, other) -> bool:
        if type(other) is not type(self):
            return False

        return (
            array_equal_shape_and_close(self.__kinetic_energy, other.__kinetic_energy)
            and array_equal_shape_and_close(
                self.__potential_energy, other.__potential_energy
            )
            and array_equal_shape_and_close(self.__total_energy, other.__total_energy)
            and array_equal_shape_and_close(self.__volume, other.__volume)
            and array_equal_shape_and_close(self.__pressure, other.__pressure)
            and array_equal_shape_and_close(self.__temperature, other.__temperature)
            and array_equal_shape_and_close(
                self.__constant_of_motion, other.__constant_of_motion
            )
            and array_equal_shape_and_close(
                self.__kinetic_energy_per_molec, other.__kinetic_energy_per_molec
            )
            and array_equal_shape_and_close(
                self.__number_of_species, other.__number_of_species
            )
        )
