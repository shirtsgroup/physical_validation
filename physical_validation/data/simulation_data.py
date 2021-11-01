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
from typing import List, Optional

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
    def compatible(data_1, data_2) -> bool:
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
    def ensemble(self) -> EnsembleData:
        r"""EnsembleData: Information on the sampled ensemble

        Returns
        -------
        ensemble : EnsembleData
        """
        return self.__ensemble

    @ensemble.setter
    def ensemble(self, ensemble: EnsembleData) -> None:
        if not isinstance(ensemble, EnsembleData):
            raise TypeError(
                "No known conversion from " + str(type(ensemble)) + "to EnsembleData"
            )
        self.__ensemble = ensemble

    @property
    def units(self) -> UnitData:
        r"""UnitsData: Information on the sampled units

        Returns
        -------
        units : UnitData
        """
        return self.__units

    @units.setter
    def units(self, units: UnitData) -> None:
        if not isinstance(units, UnitData):
            raise TypeError(
                "No known conversion from " + str(type(units)) + "to UnitData"
            )
        self.__units = units

    @property
    def observables(self) -> ObservableData:
        r"""ObservableData: Observables collected during the simulation

        Returns
        -------
        observables : ObservableData
        """
        return self.__observables

    @observables.setter
    def observables(self, observables: ObservableData) -> None:
        if not isinstance(observables, ObservableData):
            raise TypeError(
                "No known conversion from "
                + str(type(observables))
                + "to ObservableData"
            )
        self.__observables = observables

    @property
    def trajectory(self) -> TrajectoryData:
        r"""TrajectoryData: Trajectories collected during the simulation

        Returns
        -------
        trajectory : TrajectoryData
        """
        return self.__trajectory

    @trajectory.setter
    def trajectory(self, trajectory: TrajectoryData) -> None:
        if not isinstance(trajectory, TrajectoryData):
            raise TypeError(
                "No known conversion from "
                + str(type(trajectory))
                + "to TrajectoryData"
            )
        self.__trajectory = trajectory

    @property
    def system(self) -> SystemData:
        r"""SystemData: Information on the system's system

        Returns
        -------
        system : SystemData
        """
        return self.__system

    @system.setter
    def system(self, system: SystemData) -> None:
        if not isinstance(system, SystemData):
            raise TypeError(
                "No known conversion from " + str(type(system)) + "to SystemData"
            )
        self.__system = system

    @property
    def dt(self) -> float:
        r"""The timestep of the simulation run.

        Returns
        -------
        timestep : float
        """
        return self.__dt

    @dt.setter
    def dt(self, dt: float) -> None:
        dt = float(dt)
        self.__dt = dt

    def set_ensemble(
        self,
        ensemble: str,
        natoms: Optional[float] = None,
        mu: Optional[float] = None,
        volume: Optional[float] = None,
        pressure: Optional[float] = None,
        energy: Optional[float] = None,
        temperature: Optional[float] = None,
    ) -> None:
        self.__ensemble = EnsembleData(
            ensemble, natoms, mu, volume, pressure, energy, temperature
        )

    def __eq__(self, other) -> bool:
        if type(other) is not type(self):
            return False
        return (
            self.__units == other.__units
            and self.__dt == other.__dt
            and self.__system == other.__system
            and self.__ensemble == other.__ensemble
            and self.__observables == other.__observables
            and self.__trajectory == other.__trajectory
        )

    def raise_if_units_are_none(
        self,
        test_name: str,
        argument_name: str,
    ) -> None:
        r"""
        Raise if the unit data was not set

        Parameters
        ----------
        test_name
            String naming the test used for error output
        argument_name
            String naming the SimulationData argument used for error output
        """
        if self.units is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} does not contain "
                f"information about the used units.",
            )

    def raise_if_ensemble_is_invalid(
        self,
        test_name: str,
        argument_name: str,
        check_pressure: bool,
        check_mu: bool,
    ) -> None:
        r"""
        Raise InputError if the ensemble data does not hold the required
        information needed for the tests.

        Parameters
        ----------
        test_name
            String naming the test used for error output
        argument_name
            String naming the SimulationData argument used for error output
        check_pressure
            Whether to check if the pressure is defined (NPT only).
        check_mu
            Whether to check if the chemical potential is defined (muVT only).
        """
        if self.ensemble is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} does not contain "
                f"information about the sampled ensemble.",
            )
        for ensemble in "NVT", "NPT", "muVT":
            if self.ensemble.ensemble == ensemble and self.ensemble.temperature is None:
                raise pv_error.InputError(
                    argument_name,
                    f"SimulationData object provided to {test_name} indicates it was sampled "
                    f"from an {ensemble} ensemble, but does not specify the temperature.",
                )
        if (
            self.ensemble.ensemble == "NPT"
            and check_pressure
            and self.ensemble.pressure is None
        ):
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} indicates it was sampled "
                f"from an NPT ensemble, but does not specify the pressure.",
            )

        if self.ensemble.ensemble == "muVT" and check_mu and self.ensemble.mu is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} indicates it was sampled "
                f"from an muVT ensemble, but does not specify the chemical potential.",
            )

    def raise_if_system_data_is_invalid(
        self,
        test_name: str,
        argument_name: str,
        check_full_system_data_only: bool,
    ) -> None:
        r"""
        Raise InputError if the ensemble data does not hold the required
        information needed for the tests.

        Parameters
        ----------
        test_name
            String naming the test used for error output
        argument_name
            String naming the SimulationData argument used for error output
        check_full_system_data_only
            Whether to check the full system data only (number of atoms in the system,
            number of constraints in the system, number of translational and rotational
            degree of freedom reduction of the system). If false, this also checks the
            masses per atom and the molecule index and constraints.
        """
        if self.system is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} does not contain "
                f"information about the system.",
            )

        if self.system.natoms is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the number of atoms in the sampled system.",
            )
        if self.system.nconstraints is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the number of constraints in the sampled system.",
            )
        if self.system.ndof_reduction_tra is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the reduction of translational degrees of freedom in the sampled system.",
            )
        if self.system.ndof_reduction_rot is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the reduction of rotational degrees of freedom in the sampled system.",
            )

        if check_full_system_data_only:
            return

        if self.system.mass is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the mass of the atoms in the sampled system.",
            )
        if self.system.molecule_idx is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the indices of the molecules in the sampled system.",
            )
        if self.system.nconstraints_per_molecule is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} lacks information on "
                f"the number of constraints per molecule in the sampled system.",
            )

    def raise_if_observable_data_is_invalid(
        self,
        required_observables: List[str],
        test_name: str,
        argument_name: str,
    ) -> None:
        r"""
        Raise InputError if there are missing observables entries, or
        if the required observables don't have the same number of frames

        Parameters
        ----------
        required_observables
            List of strings denoting the observables required
        test_name
            String naming the test used for error output
        argument_name
            String naming the SimulationData argument used for error output
        """
        if self.observables is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} does not contain any observable data.",
            )
        for entry in required_observables:
            if self.observables[entry] is None or self.observables[entry].shape[0] == 0:
                raise pv_error.InputError(
                    argument_name,
                    f"SimulationData object provided to {test_name} lacks the observable {entry}.",
                )

        if len(required_observables) > 1:
            size_of_first_entry = self.observables[required_observables[0]].shape[0]
            if any(
                [
                    self.observables[entry].shape[0] != size_of_first_entry
                    for entry in required_observables
                ]
            ):
                raise pv_error.InputError(
                    argument_name,
                    f"SimulationData object provided to {test_name} does not have "
                    f"equal number of frames for entries "
                    + ", ".join(required_observables),
                )

    def raise_if_trajectory_data_is_invalid(
        self,
        test_name: str,
        argument_name: str,
    ) -> None:
        r"""
        Raise InputError if there are missing observables entries, or
        if the required observables don't have the same number of frames

        Parameters
        ----------
        test_name
            String naming the test used for error output
        argument_name
            String naming the SimulationData argument used for error output
        """
        if self.trajectory is None:
            raise pv_error.InputError(
                argument_name,
                f"SimulationData object provided to {test_name} does not contain any trajectory data.",
            )
        for entry in ["position", "velocity"]:
            if self.trajectory[entry] is None or self.trajectory[entry].size == 0:
                raise pv_error.InputError(
                    argument_name,
                    f"SimulationData object provided to {test_name} lacks {entry} trajectory data.",
                )
