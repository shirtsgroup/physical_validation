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
GROMOS python interface.

.. warning:: This functions in this file are not supported by the
   GROMOS development team. Any responsability is declined.
"""

import re
import gzip
import numpy as np
import warnings


class GromosInterface(object):
    @classmethod
    def read_system_from_top(cls, top_file):
        top = cls.get_blocks(top_file, unique=True)
        molecules = []

        solute_resnameblock = top['RESNAME'].split()

        solute_atomblock = top['SOLUTEATOM'].split()
        solute_natoms = int(solute_atomblock[0])
        solute_masses = []
        solute_resnames = []
        start = 1
        for natm in range(solute_natoms):
            # atnm = solute_atomblock[start + 0] - 1
            mres = solute_atomblock[start + 1]
            # panm = solute_atomblock[start + 2]
            # iac = solute_atomblock[start + 3]
            mass = solute_atomblock[start + 4]
            # cg = solute_atomblock[start + 5]
            # cgc = solute_atomblock[start + 6]
            ine = solute_atomblock[start + 7]
            ine14 = solute_atomblock[start + 8]

            start += 9 + int(ine) + int(ine14)
            solute_masses.append(float(mass))
            solute_resnames.append(solute_resnameblock[mres])

        solute_bondblock = top['BOND'].split('\n')
        solute_bondblock.pop(0)
        solute_bondhblock = top['BONDH'].split('\n')
        solute_bondhblock.pop(0)
        solute_bonds = []
        solute_bondgraph = [[] for _ in range(solute_natoms)]
        for bond in solute_bondblock.extend(solute_bondhblock):
            bond = bond.split()
            a1 = bond[0] - 1
            a2 = bond[1] - 1
            if a1 < a1:
                solute_bonds.append([a1, a2])
            else:
                solute_bonds.append([a2, a1])
            solute_bondgraph[a1].append(a2)
            solute_bondgraph[a2].append(a1)

        solute_bonds.sort(key=lambda bb: bb[0])
        first_atm = last_atm = curr_atm = 0
        curr_bonds = []
        for bond in solute_bonds:
            if bond[0] > curr_atm:
                if bond[0] > last_atm:
                    molecules.append({
                        'name': '-'.join(solute_resnames[first_atm:last_atm+1]),
                        'nmolecs': 1,
                        'natoms': first_atm - last_atm + 1,
                        'mass': solute_masses[first_atm:last_atm+1],
                        'nbonds': len(curr_bonds),
                        'bonds': curr_bonds,
                        'solvent': False
                    })
                    curr_atm = last_atm + 1
                    while curr_atm < bond[0]:
                        molecules.append({
                            'name': solute_resnames[curr_atm],
                            'nmolecs': 1,
                            'natoms': 1,
                            'mass': [solute_masses[curr_atm]],
                            'nbonds': 0,
                            'bonds': [],
                            'solvent': False
                        })
                        curr_atm += 1
                    first_atm = last_atm = curr_atm
                    curr_bonds = []
                else:
                    curr_atm = bond[0]
            if bond[1] > last_atm:
                last_atm = bond[1]
            curr_bonds.append(bond)
        else:
            molecules.append({
                'name': '-'.join(solute_resnames[first_atm:last_atm+1]),
                'nmolecs': 1,
                'natoms': first_atm - last_atm + 1,
                'mass': solute_masses[first_atm:last_atm+1],
                'nbonds': len(curr_bonds),
                'bonds': curr_bonds,
                'solvent': False
            })
            curr_atm = last_atm + 1
            while curr_atm < solute_natoms:
                molecules.append({
                    'name': solute_resnames[curr_atm],
                    'nmolecs': 1,
                    'natoms': 1,
                    'mass': [solute_masses[curr_atm]],
                    'nbonds': 0,
                    'bonds': [],
                    'solvent': False
                })
                curr_atm += 1

        solvent_atomblock = top['SOLVENTATOM'].split('\n')
        solvent_natoms = int(solvent_atomblock[0].strip())
        solvent_masses = []
        if solvent_natoms > 0:
            for line in solvent_atomblock[1:]:
                solvent_masses.append(float(line.split()[3]))
        solvent_bondblock = top['SOLVENTCONSTR'].split('\n')
        solvent_nbonds = int(solvent_bondblock[0].strip())
        solvent_bonds = []
        if solvent_nbonds > 0:
            for bond in solvent_bondblock[1:]:
                bond = bond.split()
                a1 = bond[0] - 1
                a2 = bond[1] - 1
                solvent_bonds.append([a1, a2])

        molecules.append({
            'name': 'SOL',
            'nmolecs': None,
            'natoms': solvent_natoms,
            'mass': solvent_masses,
            'nbonds': solvent_nbonds,
            'bonds': solvent_bonds,
            'solvent': True
        })

        return molecules

    @classmethod
    def read_coords(cls, files):
        traj = {}
        if isinstance(files, str):
            files = [files]
        for f in files:
            result = cls.get_blocks(f, unique=False)
            for header, content in result.items():
                if header not in traj:
                    traj[header] = []
                traj[header].extend(content)

        position = []
        velocity = []
        force = []
        box = []
        time = []

        def read_xyz(trj, reduced):
            if reduced:
                cut = 0
            else:
                cut = 24
            res = []
            for frame in trj:
                x = []
                for line in frame.split('\n'):
                    line = line[cut:].strip()
                    if not line:
                        continue
                    x.append([float(n) for n in line.split()])
                res.append(np.array(x))
            return res

        if 'POSITIONRED' in traj:
            position = read_xyz(traj['POSITIONRED'], True)
        elif 'POSITION' in traj:
            position = read_xyz(traj['POSITION'], True)
        if 'VELOCITYRED' in traj:
            velocity = read_xyz(traj['VELOCITYRED'], True)
        elif 'VELOCITY' in traj:
            velocity = read_xyz(traj['VELOCITY'], True)
        if 'FORCERED' in traj:
            force = read_xyz(traj['FORCERED'], True)
        elif 'FORCE' in traj:
            force = read_xyz(traj['FORCE'], True)
        if 'GENBOX' in traj:
            for f in traj['GENBOX']:
                f = f.split()
                if f[0] != 1:
                    warnings.warn('GENBOX boundary conditons = ' + f[0] +
                                  '. Reading from trajectory not supported - try reading volume from tre.')
                    box = []
                    break
                else:
                    box.append([float(n) for n in f[1:4]])
        if 'TIMESTEP' in traj:
            for f in traj['TIMESTEP']:
                f = f.split()
                time.append(f[1])

        result = {}
        for key, vector in zip(['position', 'velocity', 'force', 'box', 'time'],
                               [position, velocity, force, box, time]):
            vector = np.array(vector)
            if vector.size > 0:
                result[key] = vector
            else:
                result[key] = None
        return result

    @classmethod
    def get_quantities(cls, files, quantities):
        traj = {}
        for f in files:
            res = cls.get_blocks(f, unique=False)
            for header, content in res.items():
                if header not in traj:
                    traj[header] = []
                traj[header].extend(content)

        if 'ENEVERSION' in traj and traj['ENEVERSION']:
            eneversion = traj['ENEVERSION'][0].strip()
            if not all(x.strip() == eneversion for x in traj['ENEVERSION']):
                raise IOError('ENEVERSION missmatch: tre files have different versions.')
        else:
            raise IOError('ENEVERSION missing: tre files have no version.')

        blocks = cls.parse_tre_file(traj, eneversion)
        ene_ana_lib = cls.ene_ana_library()[eneversion]

        q_dict = {}
        for q in quantities:
            if q not in ene_ana_lib:
                q_dict[q] = None
                continue
            subblock = ene_ana_lib[q]['block']
            m = ene_ana_lib[q]['m']
            if 'n' not in ene_ana_lib[q]:
                values = blocks[subblock][:, m]
            else:
                n = ene_ana_lib[q]['n']
                values = blocks[subblock][:, m, n]
            if 'factor' in ene_ana_lib[q]:
                values *= ene_ana_lib[q]['factor']
            q_dict[q] = values

        return q_dict

    @classmethod
    def get_block_fields(cls, filename, unique=False):
        result = cls.get_blocks(filename, unique)
        for header in result:
            result[header] = result[header].split()
        return result

    @classmethod
    def get_block_lines(cls, filename, unique=False):
        result = cls.get_blocks(filename, unique)
        for header in result:
            result[header] = result[header].split('\n')
        return result

    @classmethod
    def get_blocks(cls, filename, unique=False):
        result = {}

        try:
            with gzip.open(filename, 'tr') as f:
                content = f.read()
        except OSError as e:
            if 'Not a gzipped file' in str(e):
                with open(filename, 'tr') as f:
                    content = f.read()
            else:
                raise e

        if content[-1] != '\n':
            content += '\n'
        content = re.sub(' *#.*\n', '\n', content)
        content = re.sub('\n+', '\n', content)
        blocks = content.split('END\n')

        for block in blocks:
            block = block.split('\n', maxsplit=1)
            header = block[0]
            if len(block) > 1:
                if header in result and unique:
                    raise IOError('Unexpected double block ' + header +
                                  'found in file ' + filename)
                elif header in result:
                    result[header].append(block[1])
                else:
                    result[header] = [block[1]]

        return result

    @classmethod
    def parse_tre_file(cls, traj, eneversion=None):

        if not eneversion or eneversion not in cls.ene_ana_library():
            raise IOError('ENEVERSION missmatch: tre file has unknown version '
                          + str(eneversion))

        ene_ana_lib = cls.ene_ana_library()[eneversion]
        enertrj = str(ene_ana_lib['ENERTRJ'])
        enertrj = re.sub(' *\n *', '\n', enertrj)
        blocks = {}
        for block in enertrj.split('\nblock'):
            block = block.strip()
            if not block:
                continue
            header = block.split('\n')[0].strip()
            frames = traj[header]
            for values in frames:
                frame = {}
                values = values.split()
                field = 0
                sizes = {}
                for line in block.split('\n')[1:]:
                    line = line.split()
                    if line[0] == 'subblock':
                        name = line[1]
                        try:
                            m = int(line[2])
                        except ValueError as e:
                            if line[2] in sizes:
                                m = sizes[line[2]]
                            else:
                                raise e
                        n = int(line[3])
                        frame[name] = []
                        if n > 1:
                            for mm in range(m):
                                row = []
                                for nn in range(n):
                                    row.append(float(values[field]))
                                    field += 1
                                frame[name].append(row)
                        else:
                            for mm in range(m):
                                frame[name].append(float(values[field]))
                                field += 1
                    if line[0] == 'size':
                        name = line[1]
                        n = int(float(values[field]))
                        sizes[name] = n
                        sizes['matrix_' + name] = int(n*(n+1)/2)
                        field += 1

                for subblock in frame:
                    if subblock in blocks:
                        blocks[subblock].append(frame[subblock])
                    else:
                        blocks[subblock] = [frame[subblock]]

        for subblock in blocks:
            blocks[subblock] = np.array(blocks[subblock])
        return blocks

    @classmethod
    def ene_ana_library(cls):
        return {
            '2015-06-23-A': {
                'kinetic_energy': {
                    'block': 'ENER',
                    'm': 1
                },
                'potential_energy': {
                    'block': 'ENER',
                    'm': 2
                },
                'total_energy': {
                    'block': 'ENER',
                    'm': 0
                },
                'volume': {
                    'block': 'VOLUME',
                    'm': 1
                },
                'pressure': {
                    'block': 'PRESSURE',
                    'm': 1,
                    'factor': 16.6057
                },
                'ENERTRJ': r"""
  block TIMESTEP
    subblock TIME 2 1
  block ENERGY03
    subblock ENER 38 1
    size NUM_BATHS
    subblock KINENER NUM_BATHS 3
    size NUM_ENERGY_GROUPS
    subblock BONDED NUM_ENERGY_GROUPS 5
    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 4
    subblock SPECIAL NUM_ENERGY_GROUPS 11
    size NUM_EDS_STATES
    subblock EDS NUM_EDS_STATES 3 
  block VOLUMEPRESSURE03
    subblock MASS 1 1
    size NUM_BATHS
    subblock TEMPERATURE  NUM_BATHS 4
    subblock VOLUME 10 1
    subblock PRESSURE 30 1
"""
            },
            '2016-01-11-A': {
                'kinetic_energy': {
                    'block': 'ENER',
                    'm': 1
                },
                'potential_energy': {
                    'block': 'ENER',
                    'm': 2
                },
                'total_energy': {
                    'block': 'ENER',
                    'm': 0
                },
                'volume': {
                    'block': 'VOLUME',
                    'm': 0
                },
                'pressure': {
                    'block': 'PRESSURE',
                    'm': 0,
                    'factor': 16.6057
                },
                'solutemp': {
                    'block': 'TEMPERATURE',
                    'm': 0,
                    'n': 0
                },
                'ENERTRJ': r"""
  block TIMESTEP
    subblock TIME 2 1
  block ENERGY03
    subblock ENER 39 1
    size NUM_BATHS
    subblock KINENER NUM_BATHS 3
    size NUM_ENERGY_GROUPS
    subblock BONDED NUM_ENERGY_GROUPS 5
    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 4
    subblock SPECIAL NUM_ENERGY_GROUPS 11
    size NUM_EDS_STATES
    subblock EDS NUM_EDS_STATES 3 
  block VOLUMEPRESSURE03
    subblock MASS 1 1
    size NUM_BATHS
    subblock TEMPERATURE  NUM_BATHS 4
    subblock VOLUME 10 1
    subblock PRESSURE 30 1
"""
            }
        }
