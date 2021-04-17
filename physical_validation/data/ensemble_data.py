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

from ..util import error as pv_error


class EnsembleData(object):
    r"""EnsembleData: Holds data defining the ensemble

    The ensemble is a string indicating the thermodynamical ensemble a simulation was
    performed in, and is any of 'NVE', 'NVT', 'NPT', 'muVT'.
    Depending on the ensemble, EnsembleData then holds additional information defining
    the ensemble, such as the number of particles N, the chemical potential mu, the
    volume V, the pressure P, the constant energy E or the temperature T. While any
    of these additional information are optional, most of them are needed by certain
    tests, such that not fully defining the ensemble results in warnings. The notable
    exception to this rule is the constant energy E for the NVE, which is not needed
    by any test and can hence be omitted without raising a warning.
    """

    @staticmethod
    def ensembles():
        return ("NVE", "NVT", "NPT", "muVT")

    def __init__(
        self,
        ensemble,
        natoms=None,
        mu=None,
        volume=None,
        pressure=None,
        energy=None,
        temperature=None,
    ):
        self.__ensemble = None
        self.__n = None
        self.__mu = None
        self.__v = None
        self.__p = None
        self.__e = None
        self.__t = None

        if ensemble not in self.ensembles():
            raise pv_error.InputError("ensemble", "Given ensemble unknown.")
        self.__ensemble = ensemble

        if ensemble == "NVE":
            if natoms is None:
                warnings.warn(ensemble + " with undefined natoms.")
            if volume is None:
                warnings.warn(ensemble + " with undefined volume.")
            # if energy is None:
            #     warnings.warn(ensemble + ' with undefined energy.')
            self.__n = natoms
            self.__v = volume
            self.__e = energy
        if ensemble == "NVT":
            if natoms is None:
                warnings.warn(ensemble + " with undefined natoms.")
            if volume is None:
                warnings.warn(ensemble + " with undefined volume.")
            if temperature is None:
                warnings.warn(ensemble + " with undefined temperature.")
            self.__n = natoms
            self.__v = volume
            self.__t = temperature
        if ensemble == "NPT":
            if natoms is None:
                warnings.warn(ensemble + " with undefined natoms.")
            if pressure is None:
                warnings.warn(ensemble + " with undefined pressure.")
            if temperature is None:
                warnings.warn(ensemble + " with undefined temperature.")
            self.__n = natoms
            self.__p = pressure
            self.__t = temperature
        if ensemble == "muVT":
            if mu is None:
                warnings.warn(ensemble + " with undefined mu.")
            if volume is None:
                warnings.warn(ensemble + " with undefined volume.")
            if temperature is None:
                warnings.warn(ensemble + " with undefined temperature.")
            self.__mu = mu
            self.__v = volume
            self.__t = temperature

    @property
    def ensemble(self):
        """Get ensemble"""
        return self.__ensemble

    @property
    def natoms(self):
        """Get natoms"""
        return self.__n

    @property
    def mu(self):
        """Get mu"""
        return self.__mu

    @property
    def volume(self):
        """Get volume"""
        return self.__v

    @property
    def pressure(self):
        """Get pressure"""
        return self.__p

    @property
    def energy(self):
        """Get energy"""
        return self.__e

    @property
    def temperature(self):
        """Get temperature"""
        return self.__t
