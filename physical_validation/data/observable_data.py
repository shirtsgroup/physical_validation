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

import numpy as np

import physical_validation.util.error as pv_error


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
        self.__nframes = -1
        self.__kinetic_energy_per_molec = None

        self.kinetic_energy = kinetic_energy
        self.potential_energy = potential_energy
        self.total_energy = total_energy
        self.volume = volume
        self.pressure = pressure
        self.temperature = temperature
        self.constant_of_motion = constant_of_motion

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

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.observables():
            raise KeyError
        return self.__getters[key](self)

    def set(self, key, value):
        self[key] = value

    def __setitem__(self, key, value):
        if key not in self.observables():
            raise KeyError
        self.__setters[key](self, value)

    @property
    def kinetic_energy(self):
        """Get kinetic_energy"""
        return self.__kinetic_energy

    @kinetic_energy.setter
    def kinetic_energy(self, kinetic_energy):
        """Set kinetic_energy"""
        if kinetic_energy is None:
            self.__kinetic_energy = None
            return
        kinetic_energy = np.array(kinetic_energy)
        if kinetic_energy.ndim != 1:
            raise pv_error.InputError("kinetic_energy", "Expected 1-dimensional array.")
        if self.nframes == -1:
            self.__nframes = kinetic_energy.size
        elif self.nframes != kinetic_energy.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__kinetic_energy = kinetic_energy

    @property
    def potential_energy(self):
        """Get potential_energy"""
        return self.__potential_energy

    @potential_energy.setter
    def potential_energy(self, potential_energy):
        """Set potential_energy"""
        if potential_energy is None:
            self.__potential_energy = None
            return
        potential_energy = np.array(potential_energy)
        if potential_energy.ndim != 1:
            raise pv_error.InputError(
                "potential_energy", "Expected 1-dimensional array."
            )
        if self.nframes == -1:
            self.__nframes = potential_energy.size
        elif self.nframes != potential_energy.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__potential_energy = potential_energy

    @property
    def total_energy(self):
        """Get total_energy"""
        return self.__total_energy

    @total_energy.setter
    def total_energy(self, total_energy):
        """Set total_energy"""
        if total_energy is None:
            self.__total_energy = None
            return
        total_energy = np.array(total_energy)
        if total_energy.ndim != 1:
            raise pv_error.InputError("total_energy", "Expected 1-dimensional array.")
        if self.nframes == -1:
            self.__nframes = total_energy.size
        elif self.nframes != total_energy.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__total_energy = total_energy

    @property
    def volume(self):
        """Get volume"""
        return self.__volume

    @volume.setter
    def volume(self, volume):
        """Set volume"""
        if volume is None:
            self.__volume = None
            return
        volume = np.array(volume)
        if volume.ndim != 1:
            raise pv_error.InputError("volume", "Expected 1-dimensional array.")
        if self.nframes == -1:
            self.__nframes = volume.size
        elif self.nframes != volume.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__volume = volume

    @property
    def pressure(self):
        """Get pressure"""
        return self.__pressure

    @pressure.setter
    def pressure(self, pressure):
        """Set pressure"""
        if pressure is None:
            self.__pressure = None
            return
        pressure = np.array(pressure)
        if pressure.ndim != 1:
            raise pv_error.InputError("pressure", "Expected 1-dimensional array.")
        if self.nframes == -1:
            self.__nframes = pressure.size
        elif self.nframes != pressure.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__pressure = pressure

    @property
    def temperature(self):
        """Get temperature"""
        return self.__temperature

    @temperature.setter
    def temperature(self, temperature):
        """Set temperature"""
        if temperature is None:
            self.__temperature = None
            return
        temperature = np.array(temperature)
        if temperature.ndim != 1:
            raise pv_error.InputError("temperature", "Expected 1-dimensional array.")
        if self.nframes == -1:
            self.__nframes = temperature.size
        elif self.nframes != temperature.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__temperature = temperature

    @property
    def constant_of_motion(self):
        """Get constant_of_motion"""
        return self.__constant_of_motion

    @constant_of_motion.setter
    def constant_of_motion(self, constant_of_motion):
        """Set constant_of_motion"""
        if constant_of_motion is None:
            self.__constant_of_motion = None
            return
        constant_of_motion = np.array(constant_of_motion)
        if constant_of_motion.ndim != 1:
            raise pv_error.InputError(
                "constant_of_motion", "Expected 1-dimensional array."
            )
        if self.nframes == -1:
            self.__nframes = constant_of_motion.size
        elif self.nframes != constant_of_motion.size:
            warnings.warn("Mismatch in number of frames. Setting `nframes = None`.")
            self.__nframes = None
        self.__constant_of_motion = constant_of_motion

    @property
    def nframes(self):
        """Get number of frames"""
        if self.__nframes is None:
            warnings.warn(
                "A mismatch in the number of frames between observables "
                "was detected. Setting `nframes = None`."
            )
        return self.__nframes

    @property
    def kinetic_energy_per_molecule(self):
        """Get kinetic_energy per molecule - used internally"""
        return self.__kinetic_energy_per_molec

    @kinetic_energy_per_molecule.setter
    def kinetic_energy_per_molecule(self, kinetic_energy):
        """Set kinetic_energy per molecule - used internally"""
        if kinetic_energy is None:
            self.__kinetic_energy_per_molec = None
            return
        # used internally - check for consistency?
        self.__kinetic_energy_per_molec = kinetic_energy

    def __eq__(self, other):
        if type(other) is not type(self):
            return False
        return (
            np.array_equal(self.__kinetic_energy, other.__kinetic_energy)
            and np.array_equal(self.__potential_energy, other.__potential_energy)
            and np.array_equal(self.__total_energy, other.__total_energy)
            and np.array_equal(self.__volume, other.__volume)
            and np.array_equal(self.__pressure, other.__pressure)
            and np.array_equal(self.__temperature, other.__temperature)
            and np.array_equal(self.__constant_of_motion, other.__constant_of_motion)
            and np.array_equal(
                self.__kinetic_energy_per_molec, other.__kinetic_energy_per_molec
            )
            and self.__nframes == other.__nframes
        )
