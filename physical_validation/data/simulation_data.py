###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (c) 2017-2021 University of Colorado Boulder               #
#              (c) 2012      The University of Virginia                   #
#                                                                         #
###########################################################################
r"""
Data structures carrying simulation data.
"""
from ..util import error as pv_error
from . import EnsembleData, ObservableData, SystemData, TrajectoryData, UnitData


class SimulationData(object):
    r"""SimulationData: System information and simulation results

    The SimulationData class holds both the information on the system and
    the results of a simulation run of that system. SimulationData contains
    all information on a simulation run needed by the physical validation
    tests. SimulationData objects can either be created directly by calling
    the class constructor, or by using a parser returning a SimulationData
    object.
    """

    @staticmethod
    def compatible(data_1, data_2):
        r"""Checks whether two simulations are compatible for common validation.

        Parameters
        ----------
        data_1 : SimulationData
        data_2 : SimulationData

        Returns
        -------
        result : bool

        """
        if not isinstance(data_1, SimulationData):
            raise pv_error.InputError("data_1", "Expected type SimulationData")
        if not isinstance(data_2, SimulationData):
            raise pv_error.InputError("data_2", "Expected type SimulationData")

        return data_1.units == data_2.units

    def __init__(
        self,
        units=None,
        dt=None,
        system=None,
        ensemble=None,
        observables=None,
        trajectory=None,
    ):
        self.__units = None
        if units is not None:
            self.units = units
        self.__dt = 0
        if dt is not None:
            self.dt = dt
        self.__system = None
        if system is not None:
            self.system = system
        self.__ensemble = None
        if ensemble is not None:
            self.ensemble = ensemble
        self.__observables = None
        if observables is not None:
            self.observables = observables
        self.__trajectory = None
        if trajectory is not None:
            self.trajectory = trajectory

    @property
    def ensemble(self):
        r"""EnsembleData: Information on the sampled ensemble

        Returns
        -------
        ensemble : EnsembleData
        """
        return self.__ensemble

    @ensemble.setter
    def ensemble(self, ensemble):
        if not isinstance(ensemble, EnsembleData):
            raise TypeError(
                "No known conversion from " + str(type(ensemble)) + "to EnsembleData"
            )
        self.__ensemble = ensemble

    @property
    def units(self):
        r"""UnitsData: Information on the sampled units

        Returns
        -------
        units : UnitData
        """
        return self.__units

    @units.setter
    def units(self, units):
        if not isinstance(units, UnitData):
            raise TypeError(
                "No known conversion from " + str(type(units)) + "to UnitData"
            )
        self.__units = units

    @property
    def observables(self):
        r"""ObservableData: Observables collected during the simulation

        Returns
        -------
        observables : ObservableData
        """
        return self.__observables

    @observables.setter
    def observables(self, observables):
        if not isinstance(observables, ObservableData):
            raise TypeError(
                "No known conversion from "
                + str(type(observables))
                + "to ObservableData"
            )
        self.__observables = observables

    @property
    def trajectory(self):
        r"""TrajectoryData: Trajectories collected during the simulation

        Returns
        -------
        trajectory : TrajectoryData
        """
        return self.__trajectory

    @trajectory.setter
    def trajectory(self, trajectory):
        if not isinstance(trajectory, TrajectoryData):
            raise TypeError(
                "No known conversion from "
                + str(type(trajectory))
                + "to TrajectoryData"
            )
        self.__trajectory = trajectory

    @property
    def system(self):
        r"""SystemData: Information on the system's system

        Returns
        -------
        system : SystemData
        """
        return self.__system

    @system.setter
    def system(self, system):
        if not isinstance(system, SystemData):
            raise TypeError(
                "No known conversion from " + str(type(system)) + "to SystemData"
            )
        self.__system = system

    @property
    def dt(self):
        r"""The timestep of the simulation run.

        Returns
        -------
        timestep : float
        """
        return self.__dt

    @dt.setter
    def dt(self, dt):
        dt = float(dt)
        self.__dt = dt

    def set_ensemble(
        self,
        ensemble,
        natoms=None,
        mu=None,
        volume=None,
        pressure=None,
        energy=None,
        temperature=None,
    ):
        self.__ensemble = EnsembleData(
            ensemble, natoms, mu, volume, pressure, energy, temperature
        )
