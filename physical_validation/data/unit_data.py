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
from typing import Optional


class UnitData(object):
    r"""UnitData: Information about the units used

    The information about units consists of different parts:

    * The name of the units (energy_str, length_str, volume_str,
                             temperature_str, pressure_str, time_str),
    * the value of kB in the used energy units, and
    * the conversion factor to GROMACS units (kJ/mol, nm, nm^3, K, bar, ps).

    The names are only used for output (console printing and plotting), and are optional.
    The conversion factors and kB are, on the other hand, used in computations and need
    to be given.
    """

    def __init__(
        self,
        kb: float,
        energy_conversion: float,
        length_conversion: float,
        volume_conversion: float,
        temperature_conversion: float,
        pressure_conversion: float,
        time_conversion: float,
        energy_str: str = "ENE",
        length_str: str = "LEN",
        volume_str: str = "VOL",
        temperature_str: str = "TEMP",
        pressure_str: str = "PRESS",
        time_str: str = "TIME",
    ):

        self.__kb = float(kb)
        self.__energy_str = str(energy_str)
        self.__energy_conversion = float(energy_conversion)
        self.__length_str = str(length_str)
        self.__length_conversion = float(length_conversion)
        self.__volume_str = str(volume_str)
        self.__volume_conversion = float(volume_conversion)
        self.__temperature_str = str(temperature_str)
        self.__temperature_conversion = float(temperature_conversion)
        self.__pressure_str = str(pressure_str)
        self.__pressure_conversion = float(pressure_conversion)
        self.__time_str = str(time_str)
        self.__time_conversion = float(time_conversion)

    @staticmethod
    def __units():
        from . import GromacsParser, LammpsParser

        return {
            "GROMACS": GromacsParser.units(),
            "LAMMPS real": LammpsParser.units("real"),
        }

    @classmethod
    def units(cls, name: Optional[str] = None):
        if name is None:
            return cls.__units().keys()

        if name in cls.__units():
            return cls.__units()[name]
        else:
            raise KeyError("Name " + name + " does not match a registered unit type.")

    def __eq__(self, other) -> bool:
        if type(other) is not type(self):
            return False

        return (
            self.kb == other.kb
            and self.energy_conversion == other.energy_conversion
            and self.length_conversion == other.length_conversion
            and self.volume_conversion == other.volume_conversion
            and self.temperature_conversion == other.temperature_conversion
            and self.pressure_conversion == other.pressure_conversion
            and self.time_conversion == other.time_conversion
        )

    @property
    def kb(self) -> float:
        """float: The value of the Boltzmann constant"""
        return self.__kb

    @property
    def energy_str(self) -> str:
        """str: Energy unit"""
        return self.__energy_str

    @property
    def length_str(self) -> str:
        """str: Length unit"""
        return self.__length_str

    @property
    def volume_str(self) -> str:
        """str: Volume unit"""
        return self.__volume_str

    @property
    def temperature_str(self) -> str:
        """str: Temperature unit"""
        return self.__temperature_str

    @property
    def pressure_str(self) -> str:
        """str: Pressure unit"""
        return self.__pressure_str

    @property
    def time_str(self) -> str:
        """str: Time unit"""
        return self.__time_str

    @property
    def energy_conversion(self) -> float:
        """float: Energy conversion factor, 1 energy_unit = energy_conversion * kJ/mol"""
        return self.__energy_conversion

    @property
    def length_conversion(self) -> float:
        """float: Length conversion factor, 1 length_unit = length_conversion * nm"""
        return self.__length_conversion

    @property
    def volume_conversion(self) -> float:
        """float: Volume conversion factor, 1 volume_unit = volume_conversion * nm^3"""
        return self.__volume_conversion

    @property
    def temperature_conversion(self) -> float:
        """float: Temperature conversion factor, 1 temperature_unit = temperature_conversion * K"""
        return self.__temperature_conversion

    @property
    def pressure_conversion(self) -> float:
        """float: Pressure conversion factor, 1 pressure_unit = pressure_conversion * bar"""
        return self.__pressure_conversion

    @property
    def time_conversion(self) -> float:
        """float: Time conversion factor, 1 time_unit = time_conversion * ps"""
        return self.__time_conversion
