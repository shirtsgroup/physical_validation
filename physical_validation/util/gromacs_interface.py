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
GROMACS python interface.

.. warning:: This is a mere place holder, as an official python API is
   currently being developed by the gromacs development team. It is
   probably neither especially elegant nor especially safe. Use of this
   module in any remotely critical application is strongly discouraged.
"""
import errno
import os
import re
import subprocess
import sys
import warnings
from typing import IO, Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..data.trajectory_data import RectangularBox
from ..util import error as pv_error

Pipe = Union[None, int, IO[Any]]


class GromacsInterface(object):
    def __init__(
        self,
        exe: str = None,
        dp: bool = False,
        includepath: Optional[Union[str, List[str]]] = None,
    ):

        self._exe = None
        self._dp = False
        self._includepath = None

        self.double = dp

        if exe is None:
            # check whether 'gmx' / 'gmx_d' is in the path
            if self._check_exe(quiet=True, exe="gmx"):
                self.exe = "gmx"
            elif self._check_exe(quiet=True, exe="gmx_d"):
                self.exe = "gmx_d"
            else:
                pv_error.InputError(
                    "exe",
                    '"gmx" and "gmx_d" not found in the path. '
                    "Set `exe` to point to a GROMACS executable.",
                )
        else:
            self.exe = exe

        if includepath is not None:
            self.includepath = includepath

    @property
    def exe(self) -> str:
        """exe is a string pointing to the gmx executable."""
        return self._exe

    @exe.setter
    def exe(self, exe: str) -> None:
        # expand '~/bin' to '/home/user/bin'
        exe = os.path.expanduser(exe)
        if self._check_exe(exe=exe):
            if os.path.dirname(exe):
                exe = os.path.abspath(exe)
            self._exe = exe

    @property
    def double(self) -> bool:
        """double is a bool defining whether the simulation was ran at double precision"""
        return self._dp

    @double.setter
    def double(self, dp: bool) -> None:
        assert isinstance(dp, bool)
        self._dp = dp

    @property
    def includepath(self) -> List[str]:
        """includepath defines a path the parser looks for system files"""
        return self._includepath

    @includepath.setter
    def includepath(self, path: Union[str, List[str]]) -> None:
        if isinstance(path, str):
            path = [path]
        self._includepath = path

    def get_quantities(
        self,
        edr: str,
        quantities: List[str],
        cwd: Optional[str] = None,
        begin: Optional[float] = None,
        end: Optional[float] = None,
        args: Optional[List[str]] = None,
    ) -> Dict[str, np.ndarray]:

        if args is None:
            args = []

        tmp_xvg = "gmxpy_" + os.path.basename(edr).replace(".edr", "") + ".xvg"
        if cwd is not None:
            tmp_xvg = os.path.join(cwd, tmp_xvg)

        q_dict = {}

        for q in quantities:
            not_found = self._create_xvg(
                edr, tmp_xvg, [q], cwd=cwd, begin=begin, end=end, args=args
            )[1]
            if q in not_found:
                q_dict[q] = None
                continue

            skip_line = re.compile("^[#,@]")
            values = []
            times = []
            with open(tmp_xvg, "r") as xvg:
                for line in xvg:
                    if skip_line.match(line):
                        continue
                    times.append(float(line.split()[0]))
                    values.append(float(line.split()[1]))

            if "time" in q_dict:
                if not np.array_equal(np.array(times), q_dict["time"]):
                    warnings.warn("Time discrepancy in " + edr)
            else:
                q_dict["time"] = np.array(times)
            q_dict[q] = np.array(values)

            os.remove(tmp_xvg)

        return q_dict

    def read_trr(
        self, trr: str
    ) -> Dict[str, Optional[Union[np.ndarray, RectangularBox]]]:
        tmp_dump = "gmxpy_" + os.path.basename(trr).replace(".trr", "") + ".dump"
        with open(tmp_dump, "w") as dump_file:
            proc = self._run(
                "dump", ["-f", trr], stdout=dump_file, stderr=subprocess.PIPE
            )
            proc.wait()

        position = []
        velocity = []
        force = []
        box = []
        with open(tmp_dump) as dump:
            x = []
            v = []
            f = []
            b = []
            for line in dump:
                if "frame" in line:
                    # new frame
                    if len(x) > 0:
                        # not the first frame - nothing to save there
                        position.append(np.array(x))
                        velocity.append(np.array(v))
                        force.append(np.array(f))
                        box.append(np.array(b))
                    x = []
                    v = []
                    f = []
                    b = []
                    continue
                if "box[" in line:
                    b.append(
                        [
                            float(value.strip())
                            for value in line.split("{", 1)[1].split("}")[0].split(",")
                        ]
                    )
                elif "x[" in line:
                    x.append(
                        [
                            float(value.strip())
                            for value in line.split("{", 1)[1].split("}")[0].split(",")
                        ]
                    )
                elif "v[" in line:
                    v.append(
                        [
                            float(value.strip())
                            for value in line.split("{", 1)[1].split("}")[0].split(",")
                        ]
                    )
                elif "f[" in line:
                    f.append(
                        [
                            float(value.strip())
                            for value in line.split("{", 1)[1].split("}")[0].split(",")
                        ]
                    )
            # end loop over file - save last arrays
            position.append(np.array(x))
            velocity.append(np.array(v))
            force.append(np.array(f))
            box.append(np.array(b))

        os.remove(tmp_dump)

        result = {}
        for key, vector in zip(
            ["position", "velocity", "force", "box"], [position, velocity, force, box]
        ):
            vector = np.array(vector)
            if vector.size > 0:
                result[key] = vector
            else:
                result[key] = None
        return result

    @staticmethod
    def read_gro(gro: str) -> Dict[str, Optional[Union[np.ndarray, RectangularBox]]]:
        with open(gro) as conf:
            x = []
            v = []
            b = []
            title = conf.readline()
            while title:
                natoms = int(conf.readline().strip())
                for _ in range(natoms):
                    line = conf.readline()[20:]
                    line = line.split()
                    x.append([float(xx) for xx in line[0:3]])
                    v.append([float(vv) for vv in line[3:6]])

                line = conf.readline()
                line = line.split()
                b.append([float(vv) for vv in line])
                title = conf.readline()

        result = {}
        for key, vector in zip(["position", "velocity", "force", "box"], [x, v, [], b]):
            vector = np.array(vector)
            if vector.size > 0:
                result[key] = vector
            else:
                result[key] = None
        return result

    @staticmethod
    def read_mdp(mdp: str) -> Dict[str, str]:
        result = {}
        with open(mdp) as f:
            for line in f:
                line = line.split(";")[0].strip()
                if not line:
                    continue
                line = line.split("=")
                # unify mdp options - all lower case, only dashes
                option = line[0].strip().replace("_", "-").lower()
                if option not in ["include", "define"]:
                    value = line[1].strip().replace("_", "-").lower()
                else:
                    value = line[1].strip()
                result[option] = value
        return result

    @staticmethod
    def write_mdp(options: Dict[str, str], mdp: str):
        with open(mdp, "w") as f:
            for key, value in options.items():
                f.write("{:24s} = {:s}\n".format(key, value))

    def read_system_from_top(
        self, top: str, define: Optional[str] = None, include: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        if not define:
            define = []
        else:
            define = [d.strip() for d in define.split("-D") if d.strip()]
        if not include:
            include = [os.getcwd(), os.path.dirname(top)]
        else:
            include = [os.getcwd(), os.path.dirname(top)] + [
                i.strip() for i in include.split("-I") if i.strip()
            ]
        superblock = None
        block = None
        nmoleculetypes = 0
        topology = {}
        with open(top) as f:
            content = self._read_top(f, include=include, define=define)

        for line in content:
            if line[0] == "[" and line[-1] == "]":
                block = line.strip("[").strip("]").strip()
                if block == "defaults" or block == "system":
                    superblock = block
                    topology[superblock] = {}
                if block == "moleculetype" or block == "molecule_type":
                    nmoleculetypes += 1
                    superblock = block + "_" + str(nmoleculetypes)
                    topology[superblock] = {}
                continue
            if superblock is None or block is None:
                raise IOError("Not a valid .top file.")
            if block in topology[superblock]:
                topology[superblock][block].append(line)
            else:
                topology[superblock][block] = [line]

        for n in range(1, nmoleculetypes + 1):
            superblock = "moleculetype_" + str(n)
            molecule = topology[superblock]["moleculetype"][0].split()[0]
            topology[molecule] = topology.pop(superblock)

        atomtype_list = topology["defaults"]["atomtypes"]
        topology["defaults"]["atomtypes"] = {}
        for atomtype in atomtype_list:
            code = atomtype.split()[0]
            topology["defaults"]["atomtypes"][code] = atomtype

        molecules = []
        for line in topology["system"]["molecules"]:
            molecule = line.split()[0]
            nmolecs = int(line.split()[1])
            natoms = len(topology[molecule]["atoms"])

            masses = []
            for atom in topology[molecule]["atoms"]:
                if len(atom.split()) >= 8:
                    masses.append(float(atom.split()[7]))
                else:
                    code = atom.split()[1]
                    masses.append(
                        float(topology["defaults"]["atomtypes"][code].split()[2])
                    )

            nbonds = 0
            nbondsh = 0
            bonds = []
            bondsh = []
            if "bonds" in topology[molecule]:
                for bond in topology[molecule]["bonds"]:
                    bond = bond.split()
                    a1 = int(bond[0]) - 1
                    a2 = int(bond[1]) - 1
                    m1 = masses[a1]
                    m2 = masses[a2]
                    if m1 > 1.008 and m2 > 1.008:
                        nbonds += 1
                        bonds.append([a1, a2])
                    else:
                        nbondsh += 1
                        bondsh.append([a1, a2])

            nangles = 0
            nanglesh = 0
            angles = []
            anglesh = []
            if "angles" in topology[molecule]:
                for angle in topology[molecule]["angles"]:
                    angle = angle.split()
                    a1 = int(angle[0]) - 1
                    a2 = int(angle[1]) - 1
                    a3 = int(angle[2]) - 1
                    m1 = masses[a1]
                    m2 = masses[a2]
                    m3 = masses[a3]
                    if m1 > 1.008 and m2 > 1.008 and m3 > 1.008:
                        nangles += 1
                        angles.append([a1, a2, a3])
                    else:
                        nanglesh += 1
                        anglesh.append([a1, a2, a3])

            settle = False
            if "settles" in topology[molecule]:
                settle = True
                bonds = []
                bondsh = [[0, 1], [0, 2], [1, 2]]

            molecules.append(
                {
                    "name": molecule,
                    "nmolecs": nmolecs,
                    "natoms": natoms,
                    "mass": masses,
                    "nbonds": [nbonds, nbondsh],
                    "bonds": bonds,
                    "bondsh": bondsh,
                    "nangles": [nangles, nanglesh],
                    "angles": angles,
                    "anglesh": anglesh,
                    "settles": settle,
                }
            )

        return molecules

    def grompp(
        self,
        mdp: str,
        top: str,
        gro: str,
        tpr: Optional[str] = None,
        cwd: str = ".",
        args: Optional[List[str]] = None,
        stdin: Optional[Pipe] = None,
        stdout: Optional[Pipe] = None,
        stderr: Optional[Pipe] = None,
    ) -> int:
        cwd = os.path.abspath(cwd)
        assert os.path.exists(os.path.join(cwd, mdp))
        assert os.path.exists(os.path.join(cwd, top))
        assert os.path.exists(os.path.join(cwd, gro))

        if args is None:
            args = []

        if tpr is None:
            tpr = os.path.basename(mdp).replace(".mdp", "") + ".tpr"
        else:
            assert os.path.exists(os.path.join(cwd, os.path.dirname(tpr)))

        args = ["-f", mdp, "-p", top, "-c", gro, "-o", tpr] + args
        proc = self._run(
            "grompp", args, cwd=cwd, stdin=stdin, stdout=stdout, stderr=stderr
        )
        proc.wait()
        return proc.returncode

    def mdrun(
        self,
        tpr: str,
        edr: Optional[str] = None,
        deffnm: Optional[str] = None,
        cwd: str = ".",
        args: Optional[List[str]] = None,
        stdin: Optional[Pipe] = None,
        stdout: Optional[Pipe] = None,
        stderr: Optional[Pipe] = None,
        mpicmd: Optional[str] = None,
    ) -> int:
        cwd = os.path.abspath(cwd)
        tpr = os.path.join(cwd, tpr)
        assert os.path.exists(cwd)
        assert os.path.exists(tpr)

        if args is None:
            args = []

        if deffnm is None:
            deffnm = os.path.basename(tpr).replace(".tpr", "")

        args = ["-s", tpr, "-deffnm", deffnm] + args
        if edr is not None:
            args += ["-e", edr]
        proc = self._run(
            "mdrun",
            args,
            cwd=cwd,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            mpicmd=mpicmd,
        )
        proc.wait()
        return proc.returncode

    def _check_exe(self, quiet: bool = False, exe: Optional[str] = None) -> bool:
        if exe is None:
            exe = self._exe
        try:
            devnull = open(os.devnull)
            exe_out = subprocess.check_output([exe, "--version"], stderr=devnull)
        except OSError as e:
            if e.errno == errno.ENOENT:
                # file not found error.
                if not quiet:
                    print("ERROR: gmx executable not found")
                    print(exe)
                return False
            else:
                raise e
        # check that output is as expected
        return bool(re.search(br":-\) GROMACS - gmx.* \(-:", exe_out))

    def _run(
        self,
        cmd: str,
        args: Optional[List[str]] = None,
        cwd: Optional[str] = None,
        stdin: Optional[Pipe] = None,
        stdout: Optional[Pipe] = None,
        stderr: Optional[Pipe] = None,
        mpicmd: Optional[str] = None,
    ) -> subprocess.Popen:
        if self.exe is None:
            raise RuntimeError(
                "Tried to use GromacsParser before setting gmx executable."
            )
        if mpicmd:
            command = [mpicmd, self.exe, cmd]
        else:
            command = [self.exe, cmd]
        command.extend(args)
        return subprocess.Popen(
            command, cwd=cwd, stdin=stdin, stdout=stdout, stderr=stderr
        )

    def _create_xvg(
        self,
        edr: str,
        xvg: str,
        quantities: List[str],
        cwd: Optional[str] = None,
        begin: Optional[float] = None,
        end: Optional[float] = None,
        args: Optional[List[str]] = None,
    ) -> Tuple[int, List[str]]:
        assert os.path.exists(edr)
        assert os.path.exists(os.path.abspath(os.path.dirname(xvg)))

        if args is None:
            args = []

        if self._dp:
            args.append("-dp")
        if begin is not None:
            args.extend(["-b", str(begin)])
        if end is not None:
            args.extend(["-e", str(end)])

        quants = ""
        for q in quantities:
            quants += str(q) + "\n"

        args = ["-f", edr, "-o", xvg] + args
        proc = self._run(
            "energy",
            args,
            cwd=cwd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        err = proc.communicate(quants.encode())[1]

        encoding = sys.stderr.encoding
        if encoding is None:
            encoding = "UTF-8"

        not_found = []
        if "does not match anything" in err.decode(encoding):
            for q in quantities:
                if "String '" + q + "' does not match anything" in err.decode(encoding):
                    not_found.append(q)

        return proc.wait(), not_found

    def _read_top(
        self, filehandler: IO, include: List[str], define: List[str]
    ) -> List[str]:
        read = [True]
        content = []
        include_dirs = include
        if self.includepath:
            include_dirs += self.includepath
        if self.exe is not None:
            # Check for standard topology locations
            try:
                full_exe = (
                    subprocess.check_output(
                        ["which", self.exe], stderr=open(os.devnull)
                    )
                    .strip()
                    .decode()
                )
            except subprocess.CalledProcessError:
                full_exe = None
            if full_exe is not None:
                # By default, when installed, exe is bin/gmx[_d],
                # while topologies are in share/gromacs/top/
                base_directory = os.path.dirname(os.path.dirname(full_exe))
                standard_install_location = os.path.join(
                    base_directory, "share", "gromacs", "top"
                )
                # When running from the build directory, exe is usually in build_folder/bin/gmx[_d],
                # while topologies are in share/top/
                root_directory = os.path.dirname(base_directory)
                standard_build_location = os.path.join(
                    os.path.join(root_directory, "share", "top")
                )
                for location in [standard_install_location, standard_build_location]:
                    if os.path.exists(location):
                        include_dirs.append(location)

        for idx, d in enumerate(include_dirs):
            # expand '~/bin' to '/home/user/bin'
            include_dirs[idx] = os.path.expanduser(d)
        for line in filehandler:
            line = line.split(";")[0].strip()
            if not line or line[0] == "*":
                continue
            if line[0] == "#":
                if line.startswith("#ifdef"):
                    option = line.replace("#ifdef", "").strip()
                    if option in define:
                        read.append(True)
                    else:
                        read.append(False)
                elif line.startswith("#ifndef"):
                    option = line.replace("#ifndef", "").strip()
                    if option not in define:
                        read.append(True)
                    else:
                        read.append(False)
                elif line.startswith("#else"):
                    read[-1] = not read[-1]
                elif line.startswith("#endif"):
                    read.pop()
                elif line.startswith("#define") and all(read):
                    option = line.replace("#define", "").strip()
                    define.append(option)
                elif line.startswith("#include") and all(read):
                    filename = (
                        line.replace("#include", "")
                        .strip()
                        .replace('"', "")
                        .replace("'", "")
                    )
                    for idir in include_dirs:
                        try:
                            ifile = open(os.path.join(idir, filename))
                            break
                        except FileNotFoundError:
                            pass
                    else:
                        msg = (
                            "Include file in .top file not found: "
                            + line
                            + "\n"
                            + "Include directories: "
                            + str(include_dirs)
                        )
                        raise IOError(msg)
                    if ifile:
                        subcontent = self._read_top(
                            ifile, [os.path.dirname(ifile.name)] + include, define
                        )
                        content.extend(subcontent)
                elif all(read):
                    raise IOError(
                        "Unknown preprocessor directive in .top file: " + line
                    )
                continue
            # end line starts with '#'

            if all(read):
                content.append(line)

        return content
