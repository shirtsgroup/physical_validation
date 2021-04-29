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
from typing import Optional

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
    def observables():
        return [
            "kinetic_energy",
            "potential_energy",
            "total_energy",
            "volume",
            "pressure",
            "temperature",
            "constant_of_motion",
        ]

    def __init__(
        self,
        kinetic_energy=None,
        potential_energy=None,
        total_energy=None,
        volume=None,
        pressure=None,
        temperature=None,
        constant_of_motion=None,
    ):
        self.__kinetic_energy = None
        self.__potential_energy = None
        self.__total_energy = None
        self.__volume = None
        self.__pressure = None
        self.__temperature = None
        self.__constant_of_motion = None
        self.__kinetic_energy_per_molec = None

        self.__getters = {
            "kinetic_energy": ObservableData.kinetic_energy.__get__,
            "potential_energy": ObservableData.potential_energy.__get__,
            "total_energy": ObservableData.total_energy.__get__,
            "volume": ObservableData.volume.__get__,
            "pressure": ObservableData.pressure.__get__,
            "temperature": ObservableData.temperature.__get__,
            "constant_of_motion": ObservableData.constant_of_motion.__get__,
        }

        self.__setters = {
            "kinetic_energy": ObservableData.kinetic_energy.__set__,
            "potential_energy": ObservableData.potential_energy.__set__,
            "total_energy": ObservableData.total_energy.__set__,
            "volume": ObservableData.volume.__set__,
            "pressure": ObservableData.pressure.__set__,
            "temperature": ObservableData.temperature.__set__,
            "constant_of_motion": ObservableData.constant_of_motion.__set__,
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

    def __getitem__(self, key):
        if key not in self.observables():
            raise KeyError
        return self.__getters[key](self)

    def __setitem__(self, key, value):
        if key not in self.observables():
            raise KeyError
        self.__setters[key](self, value)

    def __check_value(self, value, key: str) -> Optional[np.ndarray]:
        if value is None:
            return None
        value = np.array(value)
        if value.ndim != 1:
            raise pv_error.InputError(key, "Expected 1-dimensional array.")
        if self.nframes is not None and self.nframes != value.size:
            warnings.warn(
                "ObservableData: Mismatch in number of frames. Setting `nframes = None`."
            )
        return value

    @property
    def kinetic_energy(self):
        """Get kinetic_energy"""
        return self.__kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, kinetic_energy):
        """Set kinetic_energy"""
        self.__kinetic_energy = self.__check_value(kinetic_energy, "kinetic_energy")

    @property
    def potential_energy(self):
        """Get potential_energy"""
        return self.__potential_energy

    @potential_energy.setter
    def potential_energy(self, potential_energy):
        """Set potential_energy"""
        self.__potential_energy = self.__check_value(
            potential_energy, "potential_energy"
        )

    @property
    def total_energy(self):
        """Get total_energy"""
        return self.__total_energy

    @total_energy.setter
    def total_energy(self, total_energy):
        """Set total_energy"""
        self.__total_energy = self.__check_value(total_energy, "total_energy")

    @property
    def volume(self):
        """Get volume"""
        return self.__volume

    @volume.setter
    def volume(self, volume):
        """Set volume"""
        self.__volume = self.__check_value(volume, "volume")

    @property
    def pressure(self):
        """Get pressure"""
        return self.__pressure

    @pressure.setter
    def pressure(self, pressure):
        """Set pressure"""
        self.__pressure = self.__check_value(pressure, "pressure")

    @property
    def temperature(self):
        """Get temperature"""
        return self.__temperature

    @temperature.setter
    def temperature(self, temperature):
        """Set temperature"""
        self.__temperature = self.__check_value(temperature, "temperature")

    @property
    def constant_of_motion(self):
        """Get constant_of_motion"""
        return self.__constant_of_motion

    @constant_of_motion.setter
    def constant_of_motion(self, constant_of_motion):
        """Set constant_of_motion"""
        self.__constant_of_motion = self.__check_value(
            constant_of_motion, "constant_of_motion"
        )

    @property
    def nframes(self):
        """Get number of frames"""
        frames = None
        for observable in ObservableData.observables():
            if self[observable] is not None:
                if frames is not None:
                    if self[observable].size == frames:
                        continue
                    else:
                        return None
                else:
                    frames = self[observable].size

        return frames

    @property
    def kinetic_energy_per_molecule(self):
        """Get kinetic_energy per molecule - used internally"""
        return self.__kinetic_energy_per_molec

    @kinetic_energy_per_molecule.setter
    def kinetic_energy_per_molecule(self, kinetic_energy):
        """Set kinetic_energy per molecule - used internally"""
        self.__kinetic_energy_per_molec = kinetic_energy

    def __eq__(self, other):
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
        )
