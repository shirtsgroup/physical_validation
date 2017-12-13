###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
lammps_parser.py
"""
import warnings
import numpy as np
import copy

from . import parser
# py2.7 compatibility
from .simulation_data import SimulationData
from .unit_data import UnitData
from .ensemble_data import EnsembleData
from .system_data import SystemData
from .observable_data import ObservableData
from .trajectory_data import TrajectoryData
# replace lines above by this when py2.7 support is dropped:
# from . import SimulationData, UnitData, EnsembleData, SystemData, ObservableData, TrajectoryData
from ..util import error as pv_error


class LammpsParser(parser.Parser):
    """
    LammpsParser
    """
    @staticmethod
    def units(unit='lj'):
        if unit == 'real':
            return UnitData(
                kb=8.314462435405199e-3/4.18400,
                energy_str='kcal/mol',
                energy_conversion=4.18400,
                length_str='A',
                length_conversion=0.1,
                volume_str='A^3',
                volume_conversion=1e-3,
                temperature_str='K',
                temperature_conversion=1,
                pressure_str='atm',
                pressure_conversion=1.01325,
                time_str='fs',
                time_conversion=1e-3)
        else:
            raise NotImplementedError('Only LAMMPS \'unit real\' is implemented.')

    @staticmethod
    def __lammps_energy_names():
        # lammps energy codes
        return dict(kinetic_energy='KinEng',
                    potential_energy='PotEng',
                    total_energy='TotEng',
                    volume='Vol',
                    pressure='Press',
                    temperature='Temp',
                    constant_of_motion=None)

    @staticmethod
    def __default_timestep(unit='lj'):
        if unit == 'real':
            # 1 fs
            return 1
        else:
            raise NotImplementedError('Only LAMMPS \'unit real\' is implemented.')

    def __init__(self):
        super(LammpsParser, self).__init__()

    @classmethod
    def get_simulation_data(cls, ensemble=None,
                            in_file=None, data_file=None,
                            log_file=None, dump_file=None):
        r"""

        Parameters
        ----------
        ensemble: EnsembleData, optional
            As the ensemble sampled by LAMMPS can depend on many different
            options, an EnsembleData object is used
        in_file: str, optional
            A string pointing to a LAMMPS input file
        data_file: str, optional
            A string pointing to a LAMMPS data file (system information)
        log_file: str, optional
            A string pointing to a LAMMPS log file (energy trajectory)
        dump_file: str, optional
            A string pointing to a LAMMPS dump file (position / velocity trajectory)

        Returns
        -------
        result: SimulationData

        """

        # data file
        data_dict = None
        if data_file is not None:
            data_dict = cls.__read_data_file(data_file)

        # input file
        input_dict = None
        if in_file is not None:
            input_dicts = cls.__read_input_file(in_file)
            if len(input_dicts) > 1:
                warnings.warn('More than on \'run\' command found in input file.\n'
                              'Only considering last one.')
            input_dict = input_dicts[-1]

        # log file
        log_dict = None
        if log_file is not None:
            log_dicts = cls.__read_log_file(log_file)
            if len(log_dicts) > 1:
                warnings.warn('More than on run found in log file.\n'
                              'Only considering last one.')
            log_dict = log_dicts[-1]

        # dump file
        dump_dict = None
        if dump_file is not None:
            dump_dict = cls.__read_dump_file(dump_file)

        # Create SimulationData object
        result = SimulationData()

        if input_dict is not None:
            result.units = cls.units(input_dict['units'])
            if 'timestep' in input_dict:
                result.dt = float(input_dict['timestep'])
            else:
                result.dt = cls.__default_timestep(input_dict['units'])

        # Ensemble must be provided
        if ensemble is not None:
            result.ensemble = ensemble

        # trajectory data from dump
        if dump_dict is not None:
            result.trajectory = TrajectoryData(
                dump_dict['position'],
                dump_dict['velocity'])

        # system data
        if data_dict is not None:
            system = SystemData()
            system.natoms = data_dict['Header']['atoms']
            masses = [float(m[0]) for m in data_dict['Masses'].values()]
            sys_mass = []
            molecule_idx = []
            molec = -1
            for atom in data_dict['Atoms']:
                sys_mass.append(masses[atom['type']-1])
                if molec != atom['molec']:
                    molec = atom['molec']
                    molecule_idx.append(atom['n']-1)
            system.mass = sys_mass
            system.molecule_idx = molecule_idx
            system.nconstraints = 0
            system.nconstraints_per_molecule = np.zeros(len(system.molecule_idx))
            system.ndof_reduction_tra = 0
            system.ndof_reduction_rot = 0
            system.bonds = []
            system.constrained_bonds = []
            if data_dict and 'Bonds' in data_dict:
                for bond in data_dict['Bonds']:
                    system.bonds.append([a-1 for a in sorted(bond['atoms'])])
                system.bonds.sort(key=lambda bb: bb[0])
            if input_dict is not None:
                # read constraints?
                # could be fix shake / rattle or fix rigid
                constraints = []
                for style in ['shake', 'rattle']:
                    if style in input_dict['fix']:
                        if input_dict['fix'][style]['group-ID'] != 'all':
                            raise NotImplementedError('Handling of fix shake / rattle with group-ID '
                                                      'other than \'all\' not implemented.')
                        constraints.append(input_dict['fix'][style]['args'])
                if constraints:
                    const_dict = {}
                    for const_list in constraints:
                        constraint = None
                        for c in const_list[3:]:
                            if c == 'mol':
                                raise NotImplementedError('Handling of fix shake / rattle keyword '
                                                          '\'mol\' not implemented.')
                            if c in ['b', 'a', 't', 'm']:
                                constraint = c
                                if constraint not in const_dict:
                                    const_dict[constraint] = []
                                continue
                            try:
                                c_type = int(c)
                            except ValueError:
                                raise pv_error.FileFormatError(in_file,
                                                               'Parsing error in fix shake / rattle.')
                            if not constraint:
                                raise pv_error.FileFormatError(in_file,
                                                               'Parsing error in fix shake / rattle.')
                            const_dict[constraint].append(c_type)

                    if const_dict and not data_dict:
                        raise pv_error.InputError(['in_file', 'data_file'],
                                                  'in_file has constraint fix, but no data_file was given.')

                    if 'b' in const_dict and 'Bonds' in data_dict:
                        # bonds
                        for bond in data_dict['Bonds']:
                            if bond['type'] in const_dict['b']:
                                system.constrained_bonds.append(sorted(bond['atoms']))
                                molec = data_dict['Atoms'][bond['atoms'][0]-1]['molec']
                                system.nconstraints_per_molecule[molec-1] += 1
                    if 't' in const_dict and 'Bonds' in data_dict:
                        # types
                        for bond in data_dict['Bonds']:
                            a1_type = data_dict['Atoms'][bond['atoms'][0]-1]['type']
                            a2_type = data_dict['Atoms'][bond['atoms'][1]-1]['type']
                            if ((a1_type in const_dict['t'] or a2_type in const_dict['t']) and
                               sorted(bond['atoms']) not in system.constrained_bonds):
                                system.constrained_bonds.append(sorted(bond['atoms']))
                                molec = data_dict['Atoms'][bond['atoms'][0]-1]['molec']
                                system.nconstraints_per_molecule[molec-1] += 1
                    if 'm' in const_dict and 'Bonds' in data_dict:
                        # masses
                        for bond in data_dict['Bonds']:
                            a1_type = data_dict['Atoms'][bond['atoms'][0]-1]['type']
                            a2_type = data_dict['Atoms'][bond['atoms'][1]-1]['type']
                            a1_mass = masses[a1_type-1]
                            a2_mass = masses[a2_type-2]
                            for m in const_dict['m']:
                                const_mass = masses[m-1]
                                if (abs(const_mass - a1_mass) <= 0.1 or
                                   abs(const_mass - a2_mass) <= 0.1):
                                    break
                            else:
                                continue
                            if sorted(bond['atoms']) not in system.constrained_bonds:
                                system.constrained_bonds.append(sorted(bond['atoms']))
                                molec = data_dict['Atoms'][bond['atoms'][0]-1]['molec']
                                system.nconstraints_per_molecule[molec-1] += 1
                    if 'a' in const_dict and 'Angles' in data_dict:
                        # angles
                        for angle in data_dict['Angles']:
                            if (angle['type'] in const_dict['a'] and
                               sorted(angle['atoms'][:2]) in system.constrained_bonds and
                               sorted(angle['atoms'][1:]) in system.constrained_bonds):
                                molec = data_dict['Atoms'][angle['atoms'][0]-1]['molec']
                                system.nconstraints_per_molecule[molec-1] += 1
                    system.nconstraints = np.sum(system.nconstraints_per_molecule)
                    system.constrained_bonds = [[bb[0]-1, bb[1]-1] for bb in system.constrained_bonds]
                    system.constrained_bonds.sort(key=lambda bb: bb[0])
                # end if constraints

                # check for rigid bodies
                for fix_style in input_dict['fix']:
                    if fix_style.startswith('rigid'):
                        raise NotImplementedError('Handling of fix rigid not implemented.')

                # center of mass constraining
                if 'recenter' in input_dict['fix']:
                    system.ndof_reduction_tra = 3

            # end if input_dict is not None
            result.system = system

        # observable data
        if log_dict is not None:
            result.observables = ObservableData()
            lammps_energy_names = cls.__lammps_energy_names()
            if ensemble is not None and ensemble.ensemble == 'NVE':
                lammps_energy_names['constant_of_motion'] = lammps_energy_names['total_energy']
            for key, lammpskey in lammps_energy_names.items():
                if lammpskey in log_dict:
                    result.observables[key] = log_dict[lammpskey]
            if lammps_energy_names['volume'] not in log_dict:
                if dump_dict is not None:
                    vol = []
                    for b in dump_dict['box']:
                        vol.append(b[0]*b[1]*b[2])
                    result.observables['volume'] = vol

        return result

    @staticmethod
    def __read_input_file(name):
        # parse input file
        input_dicts = []
        input_dict = {'fix': {}}
        with open(name) as f:
            content = []
            previous = ''
            for line in f:
                line = line.split('#')[0].strip()
                if not line:
                    if previous:
                        content.append(previous)
                        previous = ''
                    continue
                if line[-1] == '&':
                    previous += ' ' + line[:-1]
                else:
                    content.append(previous + ' ' + line)
                    previous = ''
        for line in content:
            option = line.split(maxsplit=1)[0].strip()
            value = line.split(maxsplit=1)[1].strip()
            if option == 'run':
                input_dict['run'] = value
                input_dicts.append(copy.deepcopy(input_dict))
            elif option == 'fix':
                line = line.split()
                style = line[3]
                if style not in input_dict['fix']:
                    input_dict['fix'][style] = []
                input_dict['fix'][style] = {
                    'ID': line[1],
                    'group-ID': line[2],
                    'style': style,
                    'args': line[4:]}
            elif option == 'unfix':
                del_id = line.split()[1]
                for style in input_dict['fix']:
                    input_dict['fix'][style] = [fix for fix in input_dict['fix'][style]
                                                if fix['ID'] != del_id]
            else:
                input_dict[option] = value
        return input_dicts

    @staticmethod
    def __read_data_file(name):
        # > available blocks
        blocks = ['Header',  # 0
                  'Masses',  # 1
                  'Nonbond Coeffs',
                  'Bond Coeffs',
                  'Angle Coeffs',
                  'Dihedral Coeffs',
                  'Improper Coeffs',
                  'BondBond Coeffs',
                  'BondAngle Coeffs',
                  'MiddleBondTorsion Coeffs',
                  'EndBondTorsion Coeffs',
                  'AngleTorsion Coeffs',
                  'AngleAngleTorsion Coeffs',
                  'BondBond13 Coeffs',
                  'AngleAngle Coeffs',
                  'Atoms',  # 15
                  'Velocities',  # 16
                  'Bonds',  # 17
                  'Angles',
                  'Dihedrals',
                  'Impropers']
        file_blocks = {}

        # > read file
        with open(name) as f:
            # header section must appear first in file
            block = 'Header'
            file_blocks['Header'] = []
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
        block = 'Header'
        header_single = ['atoms',
                         'bonds',
                         'angles',
                         'dihedrals',
                         'impropers',
                         'atom types',
                         'bond types',
                         'angle types',
                         'dihedral types',
                         'improper types']
        header_double = ['xlo xhi',
                         'ylo yhi',
                         'zlo zhi']
        # default values
        data_dict[block] = {hs: 0 for hs in header_single}
        data_dict[block].update({hd: [0., 0.] for hd in header_double})
        # read out
        for line in file_blocks[block]:
            if line.split(maxsplit=1)[1] in header_single:
                hs = line.split(maxsplit=1)[1]
                data_dict[block][hs] = int(line.split(maxsplit=1)[0])
            elif line.split(maxsplit=2)[2] in header_double:
                hd = line.split(maxsplit=2)[2]
                data_dict[block][hd] = [float(line.split(maxsplit=2)[0]),
                                        float(line.split(maxsplit=2)[1])]
            else:
                raise pv_error.FileFormatError(name, 'Unknown header line')

        # > handle coeffs
        # N coeff1 coeff2 ...
        for block in blocks[1:15]:
            if block not in file_blocks:
                continue
            data_dict[block] = {}
            for line in file_blocks[block]:
                line = line.split()
                data_dict[block][int(line[0])] = line[1:]

        # > handle atoms
        # n molecule-tag atom-type q x y z nx ny nz
        block = blocks[15]
        data_dict[block] = []
        for line in file_blocks[block]:
            line = line.split()
            if len(line) == 7:
                data_dict[block].append({'n': int(line[0]),
                                         'molec': int(line[1]),
                                         'type': int(line[2]),
                                         'q': float(line[3]),
                                         'x': float(line[4]),
                                         'y': float(line[5]),
                                         'z': float(line[6])})
            else:
                data_dict[block].append({'n': int(line[0]),
                                         'molec': int(line[1]),
                                         'type': int(line[2]),
                                         'q': float(line[3]),
                                         'x': float(line[4]),
                                         'y': float(line[5]),
                                         'z': float(line[6]),
                                         'nx': float(line[7]),
                                         'ny': float(line[8]),
                                         'nz': float(line[9])})

        # > handle velocities
        # N vx vy vz
        block = blocks[16]
        if block in file_blocks:
            data_dict[block] = []
            for line in file_blocks[block]:
                line = line.split()
                data_dict[block].append({'n': int(line[0]),
                                         'vx': float(line[1]),
                                         'vy': float(line[2]),
                                         'vz': float(line[3])})

        # > handle bonds etc
        # N bond-type atom-1 atom-2 ...
        for block in blocks[17:]:
            if block not in file_blocks:
                continue
            data_dict[block] = []
            for line in file_blocks[block]:
                line = line.split()
                data_dict[block].append({'n': int(line[0]),
                                         'type': int(line[1]),
                                         'atoms': [int(c) for c in line[2:]]})

        # return dictionary
        return data_dict

    @staticmethod
    def __read_log_file(name):
        # parse log file
        def start_single(line1, line2):
            if not line1.split():
                return False
            if len(line1.split()) != len(line2.split()):
                return False
            try:
                [float(nn) for nn in line2.split()]
            except ValueError:
                return False
            return True

        def end_single(line, length):
            if len(line.split()) != length:
                return True
            try:
                [float(nn) for nn in line.split()]
            except ValueError:
                return True
            return False

        def start_multi(line):
            if '---- Step' in line and '- CPU =' in line:
                return True
            return False

        def end_multi(line):
            line = line.split()
            # right length (is it actually always 9??)
            if len(line) == 0 or len(line) % 3 != 0:
                return True
            # 2nd, 5th, 8th, ... entry must be '='
            for eq in line[1::3]:
                if eq != '=':
                    return True
            # 3rd, 6th, 9th, ... entry must be numeric
            try:
                [float(nn) for nn in line[2::3]]
            except ValueError:
                return True
            return False

        ene_trajs = []
        ene_traj = {}
        with open(name) as f:
            read_single = False
            read_multi = False
            old_line = ''
            fields = []
            for new_line in f:
                if not (read_single or read_multi):
                    if start_multi(new_line):
                        read_multi = True
                        old_line = new_line
                        continue
                    if start_single(old_line, new_line):
                        read_single = True
                        fields = old_line.split()
                        values = new_line.split()
                        for field, value in zip(fields, values):
                            ene_traj[field] = [float(value)]
                        old_line = new_line
                        continue
                    old_line = new_line
                    continue
                old_line = new_line
                if read_single:
                    if end_single(new_line, len(fields)):
                        read_single = False
                        ene_trajs.append(ene_traj)
                        ene_traj = {}
                        continue
                    for field, n in zip(fields, new_line.split()):
                        ene_traj[field].append(float(n))
                    continue
                if read_multi:
                    if end_multi(new_line):
                        if start_multi(new_line):
                            continue
                        read_multi = False
                        ene_trajs.append(ene_traj)
                        ene_traj = {}
                        continue
                    for field, n in zip(new_line.split()[0::3],
                                        new_line.split()[2::3]):
                        if field not in ene_traj:
                            ene_traj[field] = []
                        ene_traj[field].append(float(n))
                    continue

        return ene_trajs

    @staticmethod
    def __read_dump_file(name):
        # parse dump file
        # the dictionary to be filled
        dump_dict = {'position': [],
                     'velocity': [],
                     'box': []}

        # helper function checking line items
        def check_item(line_str, item):
            item = 'ITEM: ' + item
            if not line_str.startswith(item):
                raise pv_error.FileFormatError(name,
                                               'dump file: was expecting ' + item)
            return line_str.replace(item, '')

        with open(name) as f:
            line = f.readline()
            while line:
                check_item(line, 'TIMESTEP')
                f.readline()

                line = f.readline()
                check_item(line, 'NUMBER OF ATOMS')
                natoms = int(f.readline())

                line = f.readline()
                line = check_item(line, 'BOX BOUNDS')
                if len(line.split()) == 3:
                    # rectangular
                    # xx yy zz, where each of them one of
                    # p = periodic, f = fixed, s = shrink wrap,
                    # or m = shrink wrapped with a minimum value
                    line = f.readline().split()
                    bx1 = float(line[1])
                    bx2 = float(line[0])
                    bx = bx2 - bx1
                    line = f.readline().split()
                    by1 = float(line[1])
                    by2 = float(line[0])
                    by = by2 - by1
                    line = f.readline().split()
                    bz1 = float(line[1])
                    bz2 = float(line[0])
                    bz = bz2 - bz1
                elif len(line.split()) == 6:
                    # orthogonal
                    # xy xz yz xx yy zz, where xy xz yz indicates
                    # 3 tilt factors will be included, and
                    # xx yy zz being each one of
                    # p = periodic, f = fixed, s = shrink wrap,
                    # or m = shrink wrapped with a minimum value
                    raise NotImplementedError('Orthogonal box reading not implemented.')
                else:
                    raise pv_error.FileFormatError(name,
                                                   'Unknown BOX BOUNDS format')

                line = f.readline()
                line = check_item(line, 'ATOMS').split()
                if 'x' in line and 'y' in line and 'z' in line and 'id' in line:
                    scaled = False
                    iid = line.index('id')
                    irx = line.index('x')
                    iry = line.index('y')
                    irz = line.index('z')
                elif 'xs' in line and 'ys' in line and 'zs' in line and 'id' in line:
                    scaled = True
                    iid = line.index('id')
                    irx = line.index('xs')
                    iry = line.index('ys')
                    irz = line.index('zs')
                else:
                    raise pv_error.FileFormatError(name,
                                                   'Unknown ATOMS format')
                has_velocities = False
                ivx = None
                ivy = None
                ivz = None
                if 'vx' in line and 'vy' in line and 'vz' in line:
                    has_velocities = True
                    ivx = line.index('vx')
                    ivy = line.index('vy')
                    ivz = line.index('vz')

                positions = np.zeros((natoms, 3))
                if has_velocities:
                    velocities = np.zeros((natoms, 3))
                for n in range(natoms):
                    line = f.readline().split()
                    idx = int(line[iid]) - 1
                    if scaled:
                        pos = [
                            bx1 + float(line[irx])*bx,
                            by1 + float(line[iry])*by,
                            bz1 + float(line[irz])*bz
                        ]
                        positions[idx] = pos
                    else:
                        positions[idx] = [float(line[idx]) for idx in [irx, iry, irz]]
                    if has_velocities:
                        velocities[idx] = [float(line[idx]) for idx in [ivx, ivy, ivz]]

                dump_dict['position'].append(positions)
                if has_velocities:
                    dump_dict['velocity'].append(velocities)
                dump_dict['box'].append([bx, by, bz])

                # end of dump loop
                line = f.readline()

        if not dump_dict['velocity']:
            dump_dict['velocity'] = None
        return dump_dict
