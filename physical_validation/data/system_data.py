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
Data structure carrying information on the simulated system.
"""
import warnings
from typing import Dict, List, Optional

import numpy as np

from ..util import error as pv_error


class SystemData(object):
    r"""SystemData: Information about the atoms and molecules in the system.

    The information stored in SystemData objects describes the atoms and molecules
    in the system as far as the physical validation tests need it.

    The system is described in terms of

    * natoms: the total number of atoms in the system
    * nconstraints: the total number of constraints in the system
    * ndof_reduction_tra: global reduction of translational degrees of freedom (e.g.
      due to constraining the center of mass of the system)
    * ndof_reduction_rot: global reduction of rotational degrees of freedom (e.g.
      due to constraining the center of mass of the system)

    The atoms are described in terms of

    * mass: a list of the mass of every atom in the system

    The molecules are described by

    * molecule_idx: a list with the indices first atoms of every molecule (this assumes
      that the atoms are sorted by molecule)
    * nconstraints_per_molecule: a list with the number of constraints in every molecule

    Only used internally:

    * ndof_per_molecule: a list with the number of degrees of freedom of every molecule

    Reserved for future use:

    * bonds
    * constrained_bonds

    *Notes:*

    * kinetic_energy.distribution() only requires information on the system
      (natoms, nconstraints, ndof_reduction_tra, ndof_reduction_rot)
    * kinetic_energy.equipartition() additionally requires information on the atoms and molecules
      (mass, molecule_idx, nconstraints_per_molecule)

    All other tests do not require and information from SystemData.

    """

    def __init__(
        self,
        natoms: Optional[int] = None,
        nconstraints: Optional[float] = None,
        ndof_reduction_tra: Optional[float] = None,
        ndof_reduction_rot: Optional[float] = None,
        mass: Optional[np.ndarray] = None,
        molecule_idx: Optional[np.ndarray] = None,
        nconstraints_per_molecule: Optional[np.ndarray] = None,
    ):
        self.__natoms = None
        self.__nconstraints = None
        self.__ndof_reduction_tra = None
        self.__ndof_reduction_rot = None
        self.__mass = None
        self.__molecule_idx = None
        self.__nconstraints_per_molecule = None
        self.__ndof_per_molecule = None
        self.__bonds = None
        self.__constrained_bonds = None

        if natoms is not None:
            self.natoms = natoms
        if nconstraints is not None:
            self.nconstraints = nconstraints
        if ndof_reduction_tra is not None:
            self.ndof_reduction_tra = ndof_reduction_tra
        if ndof_reduction_rot is not None:
            self.ndof_reduction_rot = ndof_reduction_rot
        if mass is not None:
            self.mass = mass
        if molecule_idx is not None:
            self.molecule_idx = molecule_idx
        if nconstraints_per_molecule is not None:
            self.nconstraints_per_molecule = nconstraints_per_molecule

    @property
    def natoms(self) -> int:
        """int: Number of atoms in the system"""
        return self.__natoms

    @natoms.setter
    def natoms(self, natoms: int) -> None:
        self.__natoms = int(natoms)

    @property
    def nconstraints(self) -> float:
        """float: Total number of constraints in the system

        Does not include the reduction of degrees of freedom in the absence of
        external forces.

        """
        return self.__nconstraints

    @nconstraints.setter
    def nconstraints(self, nconstraints: float) -> None:
        self.__nconstraints = float(nconstraints)

    @property
    def ndof_reduction_tra(self) -> float:
        """float: Number of translational degrees of freedom deducted
        from 3*[# of molecules]

        """
        return self.__ndof_reduction_tra

    @ndof_reduction_tra.setter
    def ndof_reduction_tra(self, ndof_reduction_tra: float) -> None:
        self.__ndof_reduction_tra = float(ndof_reduction_tra)

    @property
    def ndof_reduction_rot(self) -> float:
        """float: Number of rotational degrees of freedom deducted
        from 3*[# of molecules]

        """
        return self.__ndof_reduction_rot

    @ndof_reduction_rot.setter
    def ndof_reduction_rot(self, ndof_reduction_rot: float) -> None:
        self.__ndof_reduction_rot = float(ndof_reduction_rot)

    @property
    def mass(self) -> np.ndarray:
        """nd-array: Mass vector for the atoms

        Setter accepts array-like objects.

        """
        return self.__mass

    @mass.setter
    def mass(self, mass: np.ndarray) -> None:
        mass = np.asarray(mass)
        if mass.ndim != 1:
            raise pv_error.InputError("mass", "Expected 1-dimensional array.")
        if self.natoms is None:
            self.natoms = mass.size
        elif mass.size != self.natoms:
            raise pv_error.InputError(
                "mass", "Mass vector does not have length == natoms."
            )
        self.__mass = mass

    @property
    def molecule_idx(self) -> np.ndarray:
        """nd-array: List of index of first atom of each molecule

        Setter accepts array-like objects.

        """
        return self.__molecule_idx

    @molecule_idx.setter
    def molecule_idx(self, molecule_idx: np.ndarray) -> None:
        molecule_idx = np.asarray(molecule_idx)
        if molecule_idx.ndim != 1:
            raise pv_error.InputError("molecule_idx", "Expected 1-dimensional array.")
        if (
            self.nconstraints_per_molecule is not None
            and self.nconstraints_per_molecule.shape != molecule_idx.shape
        ):
            warnings.warn(
                "New `molecule_idx` does not have the same"
                "shape as previously set `nconstraints_per_molecule`."
                "Setting `nconstraints_per_molecule = None` to avoid"
                "errors."
            )
            self.__nconstraints_per_molecule = None
        self.__molecule_idx = molecule_idx

    @property
    def nconstraints_per_molecule(self) -> np.ndarray:
        """nd-array: List of number of constraints per molecule

        Setter accepts array-like objects.

        """
        return self.__nconstraints_per_molecule

    @nconstraints_per_molecule.setter
    def nconstraints_per_molecule(self, nconstraints_per_molecule: np.ndarray) -> None:
        nconstraints_per_molecule = np.array(nconstraints_per_molecule)
        if nconstraints_per_molecule.ndim != 1:
            raise pv_error.InputError(
                "nconstraints_per_molecule", "Expected 1-dimensional array."
            )
        if self.molecule_idx is not None:
            if nconstraints_per_molecule.shape != self.molecule_idx.shape:
                raise pv_error.InputError(
                    "nconstraints_per_molecule",
                    "Expected `nconstraints_per_molecule` to have"
                    "the same shape as `moldecule_idx`.",
                )

        self.__nconstraints_per_molecule = nconstraints_per_molecule

    @property
    def ndof_per_molecule(self) -> Optional[List[Dict[str, float]]]:
        """nd-array: List of number of degrees of freedom per molecule"""
        return self.__ndof_per_molecule

    @ndof_per_molecule.setter
    def ndof_per_molecule(
        self, ndof_per_molecule: Optional[List[Dict[str, float]]]
    ) -> None:
        # used internally - check for consistency?
        self.__ndof_per_molecule = ndof_per_molecule

    @property
    def bonds(self) -> List[List[int]]:
        """List[List[int]]: List of bonds per molecule"""
        return self.__bonds

    @bonds.setter
    def bonds(self, bonds: List[List[int]]) -> None:
        self.__bonds = bonds

    @property
    def constrained_bonds(self) -> List[List[int]]:
        """List[List[int]]: List of constrained bonds per molecule"""
        return self.__constrained_bonds

    @constrained_bonds.setter
    def constrained_bonds(self, constrained_bonds: List[List[int]]) -> None:
        self.__constrained_bonds = constrained_bonds

    def __eq__(self, other) -> bool:
        if type(other) is not type(self):
            return False
        return (
            self.__natoms == other.__natoms
            and self.__nconstraints == other.__nconstraints
            and self.__ndof_reduction_tra == other.__ndof_reduction_tra
            and self.__ndof_reduction_rot == other.__ndof_reduction_rot
            and np.array_equal(self.__mass, other.__mass)
            and np.array_equal(self.__molecule_idx, other.__molecule_idx)
            and np.array_equal(
                self.__nconstraints_per_molecule, other.__nconstraints_per_molecule
            )
            and np.array_equal(self.__ndof_per_molecule, other.__ndof_per_molecule)
            # bonds are currently unused, so don't test for equality
            # and np.array_equal(self.__bonds, other.__bonds)
            # and np.array_equal(self.__constrained_bonds, other.__constrained_bonds)
        )
