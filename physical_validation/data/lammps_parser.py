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
lammps_parser.py
"""
from typing import Dict, List, Optional, Union

import numpy as np

from ..util import error as pv_error
from . import (
    EnsembleData,
    ObservableData,
    SimulationData,
    SystemData,
    TrajectoryData,
    UnitData,
    parser,
)


class LammpsParser(parser.Parser):
    """
    LammpsParser
    """

    @staticmethod
    def units(unit_string: str) -> UnitData:
        if unit_string == "real":
            return UnitData(
                kb=8.314462435405199e-3 / 4.18400,
                energy_str="kcal/mol",
                energy_conversion=4.18400,
                length_str="A",
                length_conversion=0.1,
                volume_str="A^3",
                volume_conversion=1e-3,
                temperature_str="K",
                temperature_conversion=1,
                pressure_str="atm",
                pressure_conversion=1.01325,
                time_str="fs",
                time_conversion=1e-3,
            )
        else:
            raise NotImplementedError("Only LAMMPS 'units real' is implemented.")

    def __init__(self):
        self.__unit = "lj"
        # lammps energy codes
        self.__lammps_energy_names = {
            "kinetic_energy": "KinEng",
            "potential_energy": "PotEng",
            "total_energy": "TotEng",
            "volume": "Vol",
            "pressure": "Press",
            "temperature": "Temp",
            "constant_of_motion": "TotEng",
        }

        # BETA warning
        print(
            "###########################################################################"
        )
        print(
            "# WARNING: The LAMMPS parser is an experimental feature under current     #"
        )
        print(
            "#          development. You can help us to improve it by reporting errors #"
        )
        print(
            "#          at https://github.com/shirtsgroup/physical_validation          #"
        )
        print(
            "#          Thank you!                                                     #"
        )
        print(
            "###########################################################################"
        )

    def get_simulation_data(
        self,
        ensemble: Optional[EnsembleData] = None,
        in_file: Optional[str] = None,
        log_file: Optional[str] = None,
        data_file: Optional[str] = None,
        dump_file: Optional[str] = None,
    ) -> SimulationData:
        """
        Parameters
        ----------
        ensemble
        in_file
        log_file
        data_file
        dump_file

        Returns
        -------
        SimulationData

        """

        # input file
        input_dict = None
        if in_file is not None:
            input_dict = self.__read_input_file(in_file)

        if input_dict is not None:
            self.__unit = input_dict["units"][0]

        # data file
        data_dict = None
        if data_file is not None:
            data_dict = self.__read_data_file(data_file)

        # log file
        log_dict = None
        if log_file is not None:
            log_dict = self.__read_log_file(log_file)

        # dump file
        dump_dict = None
        if dump_file is not None:
            dump_dict = self.__read_dump_file(dump_file)

        # Create SimulationData object
        result = SimulationData()
        result.units = self.units(self.__unit)

        # Ensemble must be provided
        if ensemble is not None:
            result.ensemble = ensemble

        # trajectory data from dump
        if dump_dict is not None:
            result.trajectory = TrajectoryData(
                dump_dict["position"], dump_dict["velocity"]
            )

        # system data
        if data_dict is not None:
            system = SystemData()
            system.natoms = data_dict["Header"]["atoms"]
            masses = data_dict["Masses"]
            mass = []
            molecule_idx = []
            molec = -1
            for atom in data_dict["Atoms"]:
                mass.append(float(masses[atom["type"]][0]))
                if molec != atom["molec"]:
                    molec = atom["molec"]
                    # LAMMPS numbering is 1-based, but internally we want 0-based
                    molecule_idx.append(atom["n"] - 1)
            system.mass = mass
            system.molecule_idx = molecule_idx
            system.nconstraints = 0
            system.nconstraints_per_molecule = np.zeros(len(system.molecule_idx))
            system.ndof_reduction_tra = 0
            system.ndof_reduction_rot = 0
            if input_dict is not None:
                if "shake" in input_dict["fix"] or "rattle" in input_dict["rattle"]:
                    print(
                        "NOTE: Found `fix shake` or `fix rattle`. Reading of\n"
                        "      constraints is currently not implemented.\n"
                        "      Please set system.nconstraints manually."
                    )
                # center of mass constraining
                if "recenter" in input_dict["fix"]:
                    system.ndof_reduction_tra = 3
            result.system = system

        # observable data
        if log_dict is not None:
            result.observables = ObservableData()
            for key, lammpskey in self.__lammps_energy_names.items():
                if lammpskey in log_dict:
                    result.observables[key] = log_dict[lammpskey]
            if self.__lammps_energy_names["volume"] not in log_dict:
                if dump_dict is not None:
                    vol = []
                    for b in dump_dict["box"]:
                        vol.append(b[0] * b[1] * b[2])
                    if len(vol) == 1:
                        vol = vol * result.observables.nframes
                    if len(vol) != result.observables.nframes and np.allclose(
                        [vol[0]] * len(vol), vol
                    ):
                        vol = [vol[0]] * result.observables.nframes
                    key = "volume"
                    result.observables[key] = vol

        return result

    @staticmethod
    def __read_input_file(
        name: str,
    ) -> Dict[str, Union[str, List[str], List[Dict[str, Union[str, List[str]]]]]]:
        # parse input file
        input_dict = {}
        with open(name) as f:
            for line in f:
                line = line.split("#")[0].strip()
                if not line:
                    continue
                option = line.split(maxsplit=1)[0].strip()
                value = line.split(maxsplit=1)[1].strip()
                if option == "fix":
                    if "fix" not in input_dict:
                        input_dict["fix"] = {}
                    line = line.split()
                    style = line[3]
                    if style not in input_dict["fix"]:
                        input_dict["fix"][style] = []
                    input_dict["fix"][style].append(
                        {
                            "ID": line[1],
                            "group-ID": line[2],
                            "style": style,
                            "args": line[4:],
                        }
                    )
                elif option == "unfix":
                    del_id = line.split()[1]
                    for style in input_dict["fix"]:
                        input_dict["fix"][style] = [
                            fix
                            for fix in input_dict["fix"][style]
                            if fix["ID"] != del_id
                        ]
                elif option in input_dict:
                    input_dict[option].append(value)
                else:
                    input_dict[option] = [value]
        return input_dict

    @staticmethod
    def __read_data_file(
        name: str,
    ) -> Dict[
        str,
        Union[
            Dict[str, Union[float, List[float]]],
            Dict[int, List[Union[str, float]]],
            List[Dict[str, Union[int, float]]],
        ],
    ]:
        # > available blocks
        blocks = [
            "Header",  # 0
            "Masses",  # 1
            "Nonbond Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
            "BondBond Coeffs",
            "BondAngle Coeffs",
            "MiddleBondTorsion Coeffs",
            "EndBondTorsion Coeffs",
            "AngleTorsion Coeffs",
            "AngleAngleTorsion Coeffs",
            "BondBond13 Coeffs",
            "AngleAngle Coeffs",
            "Atoms",  # 15
            "Velocities",  # 16
            "Bonds",  # 17
            "Angles",
            "Dihedrals",
            "Impropers",
        ]
        file_blocks = {}

        # > read file
        with open(name) as f:
            # header section must appear first in file
            block = "Header"
            file_blocks["Header"] = []
            # 1st 2 lines are ignored
            next(f)
            next(f)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line in blocks:
                    block = line
                    file_blocks[block] = []
                    continue
                file_blocks[block].append(line)

        data_dict = {}

        # > handle header
        block = "Header"
        header_single = [
            "atoms",
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
            "atom types",
            "bond types",
            "angle types",
            "dihedral types",
            "improper types",
        ]
        header_double = ["xlo xhi", "ylo yhi", "zlo zhi"]
        # default values
        # Dict[str, float]
        data_dict[block] = {hs: 0 for hs in header_single}
        # Dict[str, List[float]]
        data_dict[block].update({hd: [0.0, 0.0] for hd in header_double})
        # read out
        for line in file_blocks[block]:
            if line.split(maxsplit=1)[1] in header_single:
                hs = line.split(maxsplit=1)[1]
                data_dict[block][hs] = int(line.split(maxsplit=1)[0])
            elif line.split(maxsplit=2)[2] in header_double:
                hs = line.split(maxsplit=2)[2]
                data_dict[block][hs] = [
                    float(line.split(maxsplit=2)[0]),
                    float(line.split(maxsplit=2)[1]),
                ]
            else:
                raise pv_error.FileFormatError(name, "Unknown header line")

        # > handle coeffs
        # N type coeff1 coeff2 ...
        for block in blocks[1:15]:
            if block not in file_blocks:
                continue
            data_dict[block] = {}
            for line in file_blocks[block]:
                line = line.split()
                # Dict[int, List[Union[str, float]]]
                data_dict[block][int(line[0])] = [line[1]] + [
                    float(c) for c in line[2:]
                ]

        # > handle atoms
        # n molecule-tag atom-type q x y z nx ny nz
        block = blocks[15]
        data_dict[block] = []
        for line in file_blocks[block]:
            line = line.split()
            # List[Dict[str, Union[int, float]]]
            if len(line) == 7:
                data_dict[block].append(
                    {
                        "n": int(line[0]),
                        "molec": int(line[1]),
                        "type": float(line[2]),
                        "q": float(line[3]),
                        "x": float(line[4]),
                        "y": float(line[5]),
                        "z": float(line[6]),
                    }
                )
            else:
                data_dict[block].append(
                    {
                        "n": int(line[0]),
                        "molec": int(line[1]),
                        "type": float(line[2]),
                        "q": float(line[3]),
                        "x": float(line[4]),
                        "y": float(line[5]),
                        "z": float(line[6]),
                        "nx": float(line[7]),
                        "ny": float(line[8]),
                        "nz": float(line[9]),
                    }
                )

        # > handle velocities
        # N vx vy vz
        block = blocks[16]
        if block in file_blocks:
            data_dict[block] = []
            for line in file_blocks[block]:
                line = line.split()
                data_dict[block].append(
                    {
                        "n": int(line[0]),
                        "vx": float(line[1]),
                        "vy": float(line[2]),
                        "vz": float(line[3]),
                    }
                )

        # > handle bonds etc
        # N bond-type atom-1 atom-2 ...
        for block in blocks[17:]:
            if block not in file_blocks:
                continue
            data_dict[block] = []
            for line in file_blocks[block]:
                line = line.split()
                data_dict[block].append(
                    {"n": int(line[0]), "atoms": [int(c) for c in line[1:]]}
                )

        # return dictionary
        return data_dict

    @staticmethod
    def __read_log_file(name: str) -> Dict[str, List[float]]:
        # parse log file
        def start_single(line1: str, line2: str) -> bool:
            if not line1.split():
                return False
            if len(line1.split()) != len(line2.split()):
                return False
            try:
                [float(nn) for nn in line2.split()]
            except ValueError:
                return False
            return True

        def end_single(line: str, length: int) -> bool:
            if len(line.split()) != length:
                return True
            try:
                [float(nn) for nn in line.split()]
            except ValueError:
                return True
            return False

        def start_multi(line: str) -> bool:
            if "---- Step" in line and "- CPU =" in line:
                return True
            return False

        def end_multi(line: str) -> bool:
            line = line.split()
            # right length (is it actually always 9??)
            if len(line) == 0 or len(line) % 3 != 0:
                return True
            # 2nd, 5th, 8th, ... entry must be '='
            for eq in line[1::3]:
                if eq != "=":
                    return True
            # 3rd, 6th, 9th, ... entry must be numeric
            try:
                [float(nn) for nn in line[2::3]]
            except ValueError:
                return True
            return False

        ene_traj = {}
        nreads = 0
        with open(name) as f:
            read_single = False
            read_multi = False
            continued = False
            old_line = ""
            fields = []
            for new_line in f:
                if read_single:
                    if end_single(new_line, len(fields)):
                        read_single = False
                        continued = True
                    else:
                        for field, n in zip(fields, new_line.split()):
                            ene_traj[field].append(float(n))
                if read_multi:
                    if end_multi(new_line):
                        read_multi = False
                        continued = True
                    else:
                        for field, n in zip(
                            new_line.split()[0::3], new_line.split()[2::3]
                        ):
                            if field not in ene_traj:
                                ene_traj[field] = []
                            ene_traj[field].append(float(n))

                if not (read_single or read_multi):
                    if start_multi(new_line):
                        if not continued:
                            ene_traj = {}
                            nreads += 1
                        read_multi = True
                        old_line = new_line
                    if start_single(old_line, new_line):
                        if not continued:
                            ene_traj = {}
                            nreads += 1
                        read_single = True
                        fields = new_line.split()
                        for field in fields:
                            if field not in ene_traj:
                                ene_traj[field] = []

                old_line = new_line
                continued = False
        if nreads > 1:
            print(
                "NOTE: Multiple runs found in log file. Assumed prior runs\n"
                "      were equilibration runs and used only last run."
            )

        return ene_traj

    @staticmethod
    def __read_dump_file(name: str) -> Dict[str, List[List[Union[float, List[float]]]]]:
        # parse dump file
        # the dictionary to be filled
        dump_dict = {"position": [], "velocity": [], "box": []}

        # helper function checking line items
        def check_item(line_str: str, item: str) -> str:
            item = "ITEM: " + item
            if not line_str.startswith(item):
                raise pv_error.FileFormatError(name, "dump file: was expecting " + item)
            return line_str.replace(item, "")

        with open(name) as f:
            line = f.readline()
            while line:
                check_item(line, "TIMESTEP")
                f.readline()
                line = f.readline()
                check_item(line, "NUMBER OF ATOMS")
                natoms = int(f.readline())

                line = f.readline()
                line = check_item(line, "BOX BOUNDS")
                bx = 0
                by = 0
                bz = 0
                if len(line.split()) == 3:
                    # rectangular
                    # xx yy zz, where each of them one of
                    # p = periodic, f = fixed, s = shrink wrap,
                    # or m = shrink wrapped with a minimum value
                    line = f.readline().split()
                    bx = float(line[1]) - float(line[0])
                    line = f.readline().split()
                    by = float(line[1]) - float(line[0])
                    line = f.readline().split()
                    bz = float(line[1]) - float(line[0])
                elif len(line.split()) == 6:
                    # orthogonal
                    # xy xz yz xx yy zz, where xy xz yz indicates
                    # 3 tilt factors will be included, and
                    # xx yy zz being each one of
                    # p = periodic, f = fixed, s = shrink wrap,
                    # or m = shrink wrapped with a minimum value
                    raise NotImplementedError("Orthogonal box reading not implemented.")

                line = f.readline()
                line = check_item(line, "ATOMS").split()
                if "x" not in line or "y" not in line or "z" not in line:
                    raise pv_error.FileFormatError(name, "No positions in dump file.")
                irx = line.index("x")
                iry = line.index("y")
                irz = line.index("z")
                has_velocities = False
                ivx = None
                ivy = None
                ivz = None
                if "vx" in line and "vy" in line and "vz" in line:
                    has_velocities = True
                    ivx = line.index("vx")
                    ivy = line.index("vy")
                    ivz = line.index("vz")

                positions = []
                velocities = []
                for n in range(natoms):
                    line = f.readline().split()
                    positions.append([float(line[idx]) for idx in [irx, iry, irz]])
                    if has_velocities:
                        velocities.append([float(line[idx]) for idx in [ivx, ivy, ivz]])

                dump_dict["position"].append(positions)
                dump_dict["velocity"].append(velocities)
                dump_dict["box"].append([bx, by, bz])

                # end of dump loop
                line = f.readline()

        return dump_dict
