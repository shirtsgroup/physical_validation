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
gromacs_parser.py
"""
import warnings
from typing import List, Optional, Union

import numpy as np

from ..util import error as pv_error
from ..util.gromacs_interface import GromacsInterface
from . import (
    EnsembleData,
    ObservableData,
    SimulationData,
    SystemData,
    TrajectoryData,
    UnitData,
    parser,
)
from .trajectory_data import RectangularBox


class GromacsParser(parser.Parser):
    """
    GromacsParser
    """

    @staticmethod
    def units() -> UnitData:
        # Gromacs uses kJ/mol
        return UnitData(
            kb=8.314462435405199e-3,
            energy_str="kJ/mol",
            energy_conversion=1.0,
            length_str="nm",
            length_conversion=1.0,
            volume_str="nm^3",
            volume_conversion=1.0,
            temperature_str="K",
            temperature_conversion=1.0,
            pressure_str="bar",
            pressure_conversion=1.0,
            time_str="ps",
            time_conversion=1.0,
        )

    def __init__(
        self, exe: Optional[str] = None, includepath: Union[str, List[str]] = None
    ):
        r"""
        Create a GromacsParser object

        Parameters
        ----------
        exe: str, optional
            Path to a gmx executable (or simply the executable name, if it is in the path)
            Default: Looks for `gmx`, then for `gmx_d` in the path. If neither is found, `exe` is
                     set to None, and any parsing including simulation trajectories (`edr`, `trr`
                     and `gro` arguments in `get_simulation_data()`) will fail.
        includepath: str or List[str], optional
            Path or list of paths to location(s) of topology file. Is used for the lookup of
            `#include` statements in topologies.
            Default: None - no additional topology location. Lookup will be restricted to current
                     directory and location of the `top` file given to `get_simulation_data()`,
                     plus any include locations added to the `mdp` file.
        """
        super(GromacsParser, self).__init__()
        self.__interface = GromacsInterface(exe=exe, includepath=includepath)
        # gmx energy codes
        self.__gmx_energy_names = {
            "kinetic_energy": "Kinetic-En.",
            "potential_energy": "Potential",
            "total_energy": "Total-Energy",
            "volume": "Volume",
            "pressure": "Pressure",
            "temperature": "Temperature",
            "constant_of_motion": "Conserved-En.",
        }

    def get_simulation_data(
        self,
        mdp: Optional[str] = None,
        top: Optional[str] = None,
        edr: Optional[str] = None,
        trr: Optional[str] = None,
        gro: Optional[str] = None,
    ) -> SimulationData:
        r"""

        Parameters
        ----------
        mdp: str, optional
            A string pointing to a .mdp file
        top: str, optional
            A string pointing to a .top file
        edr: str, optional
            A string pointing to a .edr file
        trr: str, optional
            A string pointing to a .trr file
        gro: str, optional
            A string pointing to a .gro file (Note: if also trr is given, gro is ignored)

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the results of the simulation as described by
            the provided GROMACS files.

        """
        result = SimulationData()
        result.units = self.units()

        # trajectories (might be used later for the box...)
        trajectory_dict = None
        if trr is not None:
            if gro is not None:
                warnings.warn("`trr` and `gro` given. Ignoring `gro`.")

            trajectory_dict = self.__interface.read_trr(trr)

            # check box shape
            if trajectory_dict["box"].ndim == 2:
                if (
                    trajectory_dict["box"] - np.diag(np.diag(trajectory_dict["box"]))
                ).any():
                    raise NotImplementedError("Triclinic boxes not implemented.")
                else:
                    box = RectangularBox(np.diag(trajectory_dict["box"]))
            elif trajectory_dict["box"].ndim == 3:
                if np.array(
                    [b - np.diag(np.diag(b)) for b in trajectory_dict["box"]]
                ).any():
                    raise NotImplementedError("Triclinic boxes not implemented.")
                else:
                    box = RectangularBox(
                        np.array([np.diag(b) for b in trajectory_dict["box"]])
                    )
            else:
                raise RuntimeError("Unknown box shape.")
            trajectory_dict["box"] = box
        elif gro is not None:
            trajectory_dict = self.__interface.read_gro(gro)

            # check box shape
            if trajectory_dict["box"].ndim == 1:
                if (
                    trajectory_dict["box"].size > 3
                    and (trajectory_dict["box"][3:]).any()
                ):
                    raise NotImplementedError("Triclinic boxes not implemented.")
                else:
                    box = RectangularBox(trajectory_dict["box"][:3])
            elif trajectory_dict["box"].ndim == 2:
                if (
                    trajectory_dict["box"].shape[1] > 3
                    and (trajectory_dict["box"][:, 3:]).any()
                ):
                    raise NotImplementedError("Triclinic boxes not implemented.")
                else:
                    box = RectangularBox(trajectory_dict["box"][:, :3])
            else:
                raise RuntimeError("Unknown box shape.")
            trajectory_dict["box"] = box

        # simulation parameters & system
        if mdp is not None and top is not None:
            mdp_options = self.__interface.read_mdp(mdp)
            define = None
            include = None
            if "define" in mdp_options:
                define = mdp_options["define"]
            if "include" in mdp_options:
                include = mdp_options["include"]
            molecules = self.__interface.read_system_from_top(
                top, define=define, include=include
            )

            if "dt" in mdp_options:
                result.dt = float(mdp_options["dt"])

            natoms = 0
            mass = []
            constraints_per_molec = []
            angles = (
                "constraints" in mdp_options
                and mdp_options["constraints"] == "all-angles"
            )
            angles_h = (
                angles
                or "constraints" in mdp_options
                and mdp_options["constraints"] == "h-angles"
            )
            bonds = (
                angles_h
                or "constraints" in mdp_options
                and mdp_options["constraints"] == "all-bonds"
            )
            bonds_h = (
                bonds
                or "constraints" in mdp_options
                and mdp_options["constraints"] == "h-bonds"
            )

            molecule_idx = []
            next_molec = 0
            molec_bonds = []
            molec_bonds_constrained = []
            for molecule in molecules:
                natoms += molecule["nmolecs"] * molecule["natoms"]
                for n in range(0, molecule["nmolecs"]):
                    molecule_idx.append(next_molec)
                    next_molec += molecule["natoms"]
                mass.extend(molecule["mass"] * molecule["nmolecs"])
                constraints = 0
                constrained_bonds = []
                all_bonds = molecule["bonds"] + molecule["bondsh"]
                if molecule["settles"]:
                    constraints = 3
                    constrained_bonds = all_bonds
                else:
                    if bonds:
                        constraints += molecule["nbonds"][0]
                        constrained_bonds.extend(molecule["bonds"])
                    if bonds_h:
                        constraints += molecule["nbonds"][1]
                        constrained_bonds.extend(molecule["bondsh"])
                    if angles:
                        constraints += molecule["nangles"][0]
                    if angles_h:
                        constraints += molecule["nangles"][1]
                constraints_per_molec.extend([constraints] * molecule["nmolecs"])
                molec_bonds.extend([all_bonds] * molecule["nmolecs"])
                molec_bonds_constrained.extend(
                    [constrained_bonds] * molecule["nmolecs"]
                )

            system = SystemData()
            system.natoms = natoms
            system.mass = mass
            system.molecule_idx = molecule_idx
            system.nconstraints = np.sum(constraints_per_molec)
            system.nconstraints_per_molecule = constraints_per_molec
            system.ndof_reduction_tra = 3
            system.ndof_reduction_rot = 0
            if "comm-mode" in mdp_options:
                if mdp_options["comm-mode"] == "linear":
                    system.ndof_reduction_tra = 3
                elif mdp_options["comm-mode"] == "angular":
                    system.ndof_reduction_tra = 3
                    system.ndof_reduction_rot = 3
                if mdp_options["comm-mode"] == "none":
                    system.ndof_reduction_tra = 0
            system.bonds = molec_bonds
            system.constrained_bonds = molec_bonds_constrained
            result.system = system

            if trajectory_dict is not None:
                # now that we know the bonds, we can gather & save the trajectory
                trajectory_dict["position"] = trajectory_dict["box"].gather(
                    trajectory_dict["position"], molec_bonds, molecule_idx
                )
                result.trajectory = TrajectoryData(
                    trajectory_dict["position"], trajectory_dict["velocity"]
                )

            thermostat = (
                "tcoupl" in mdp_options
                and mdp_options["tcoupl"]
                and mdp_options["tcoupl"] != "no"
            )
            stochastic_dyn = "integrator" in mdp_options and mdp_options[
                "integrator"
            ] in ["sd", "sd2", "bd"]
            constant_temp = thermostat or stochastic_dyn
            temperature = None
            if constant_temp:
                ref_t = [float(t) for t in mdp_options["ref-t"].split()]
                if len(ref_t) == 1 or np.allclose(ref_t, [ref_t[0]] * len(ref_t)):
                    temperature = ref_t[0]
                else:
                    raise pv_error.InputError(
                        "mdp",
                        "Ensemble definition ambiguous: Different t-ref values found.",
                    )

            constant_press = (
                "pcoupl" in mdp_options
                and mdp_options["pcoupl"]
                and mdp_options["pcoupl"] != "no"
            )
            volume = None
            pressure = None
            if constant_press:
                ref_p = [float(p) for p in mdp_options["ref-p"].split()]
                if len(ref_p) == 1 or np.allclose(ref_p, [ref_p[0]] * len(ref_p)):
                    pressure = ref_p[0]
                else:
                    raise pv_error.InputError(
                        "mdp",
                        "Ensemble definition ambiguous: Different p-ref values found.",
                    )
            else:
                if trajectory_dict is not None:
                    box = trajectory_dict["box"].box[0]
                    # Different box shapes?
                    if box.ndim == 1:
                        volume = box[0] * box[1] * box[2]
                    elif box.ndim == 2:
                        volume = box[0, 0] * box[1, 1] * box[2, 2]

            if constant_temp and constant_press:
                ens = "NPT"
            elif constant_temp:
                ens = "NVT"
            else:
                ens = "NVE"

            if ens == "NVE":
                self.__gmx_energy_names["constant_of_motion"] = "Total-Energy"
            else:
                self.__gmx_energy_names["constant_of_motion"] = "Conserved-En."

            result.ensemble = EnsembleData(
                ens,
                natoms=natoms,
                volume=volume,
                pressure=pressure,
                temperature=temperature,
            )
        elif trajectory_dict is not None:
            # we don't know the system, so we can't gather, but save it anyway
            result.trajectory = TrajectoryData(
                trajectory_dict["position"], trajectory_dict["velocity"]
            )

        if edr is not None:
            observable_dict = self.__interface.get_quantities(
                edr, list(self.__gmx_energy_names.values()), args=["-dp"]
            )

            # constant volume simulations don't write out the volume in .edr file
            if (
                observable_dict["Volume"] is None
                and result.ensemble is not None
                and result.ensemble.volume is not None
            ):
                nframes = observable_dict["Pressure"].size
                observable_dict["Volume"] = np.ones(nframes) * result.ensemble.volume

            result.observables = ObservableData()
            for key, gmxkey in self.__gmx_energy_names.items():
                result.observables[key] = observable_dict[gmxkey]

        return result
