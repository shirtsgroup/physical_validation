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
gromos_parser.py
"""
import warnings

import numpy as np

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
from ..util.gromos_interface import GromosInterface


class GromosParser(parser.Parser):
    """
    GromosParser
    """

    @staticmethod
    def units():
        # Note: for pressure, Gromos is internally using kJ/mol/nm^3, but
        # GromosInterface does the conversion to bar
        return UnitData(
            kb=8.314462435405199e-3,
            energy_str='kJ/mol',
            energy_conversion=1.0,
            length_str='nm',
            length_conversion=1.0,
            volume_str='nm^3',
            volume_conversion=1.0,
            temperature_str='K',
            temperature_conversion=1.0,
            pressure_str='bar',
            pressure_conversion=1.0,
            time_str='ps',
            time_conversion=1.0)

    def __init__(self):
        super(GromosParser, self).__init__()

    def get_simulation_data(self,
                            imd=None, top=None,
                            trc=None, trv=None,
                            cnf=None, tre=None):
        r"""

        Parameters
        ----------
        imd: str or List[str], optional
            A string or a list of strings pointing to a .imd file
        top: str or List[str], optional
            A string or a list of strings pointing to a .top file
        trc: str or List[str], optional
            A string or a list of strings pointing to a .trc file
        trv: str or List[str], optional
            A string or a list of strings pointing to a .trv file
        cnf: str or List[str], optional
            A string or a list of strings pointing to a .cnf file
            (Note: if also trc is given, cnf is ignored)
        tre: str or List[str], optional
            A string or a list of strings pointing to a .tre file

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the result of the simulation
            as provided by the input files.

        """
        result = SimulationData()
        result.units = self.units()

        trajectory_dict = {}
        if trc:
            if cnf:
                warnings.warn('`trc` and `cnf` given. Ignoring `cnf`.')
            trc_dict = GromosInterface.read_coords(trc)
            trajectory_dict['position'] = trc_dict['position']
            if trc_dict['box']:
                trajectory_dict['volume'] = np.prod(trc_dict['box'], axis=1)
            if trv:
                trv_dict = GromosInterface.read_coords(trv)
                if not np.allclose(trc_dict['time'], trv_dict['time']):
                    warnings.warn('`trc` and `trv` have different time information.')
                trajectory_dict['velocity'] = trv_dict['velocity']
            else:
                trajectory_dict['velocity'] = None
        elif cnf:
            cnf_dict = GromosInterface.read_coords(cnf)
            trajectory_dict['position'] = cnf_dict['position']
            trajectory_dict['velocity'] = cnf_dict['velocity']
            if cnf_dict['box']:
                trajectory_dict['volume'] = np.prod(cnf_dict['box'], axis=1)
        else:
            trajectory_dict['position'] = None
            trajectory_dict['velocity'] = None

        result.trajectory = TrajectoryData(
            trajectory_dict['position'],
            trajectory_dict['velocity'])

        observable_dict = {}
        if tre:
            observables = [
                'kinetic_energy',
                'potential_energy',
                'total_energy',
                'volume',
                'pressure'
            ]
            observable_dict = GromosInterface.get_quantities(tre, observables)

            result.observables = ObservableData()
            for key in observables:
                result.observables[key] = observable_dict[key]

        molecules = []
        if top:
            molecules = GromosInterface.read_system_from_top(top)

        system = SystemData()
        if imd:
            input_dict = GromosInterface.get_block_fields(imd, unique=True)
            nsolvent_molecules = int(input_dict['SYSTEM'][1])

            # check constraints
            constbonds = False
            constbondsh = False
            if 'CONSTRAINTS' in input_dict:
                ntc = int(input_dict['CONSTRAINTS'][0])
                if ntc == 4:
                    raise NotImplementedError('Handling of CONSTRAINT option "specified" not implemented.')
                if ntc == 2:
                    constbondsh = True
                if ntc == 3:
                    constbonds = True

            # count atoms & constraints
            natoms = 0
            masses = []
            molecule_idx = [0]
            nconstraints_per_molec = []
            all_bonds = []
            constrained_bonds = []
            if molecules:
                for molecule in molecules:
                    # natoms
                    if molecule['solvent']:
                        molecule['nmolecs'] = nsolvent_molecules
                    natoms += molecule['nmolecs'] * molecule['natoms']

                    # masses
                    masses.extend(molecule['masses'] * molecule['nmolecs'])

                    # molecule_idx
                    for _ in range(molecule['nmolecs']):
                        molecule_idx.append(molecule_idx[-1] + molecule['natoms'])

                    # nconstraints, bonds
                    if molecule['solvent']:
                        nconstraints_per_molec.extend([molecule['nbonds']] * molecule['nmolecs'])
                        all_bonds.extend(molecule['bonds'] * molecule['nmolecs'])
                        constrained_bonds.extend(molecule['bonds'] * molecule['nmolecs'])
                    else:
                        all_bonds.extend(molecule['bonds'] * molecule['nmolecs'])
                        all_bonds.extend(molecule['bondsh'] * molecule['nmolecs'])
                        nconstraints = 0
                        if constbondsh:
                            nconstraints += molecule['nbonds'][1]
                            constrained_bonds.extend(molecule['bondsh'] * molecule['nmolecs'])
                        if constbonds:
                            nconstraints += molecule['nbonds'][0]
                            constrained_bonds.extend(molecule['bonds'] * molecule['nmolecs'])
                        nconstraints_per_molec.extend([nconstraints] * molecule['nmolecs'])

            # end if molecules
            molecule_idx.pop()  # remove last entry (no more molecule *starting* there)

            system.natoms = natoms
            system.nconstraints = np.sum(nconstraints_per_molec)
            system.ndof_reduction_tra = 0
            system.ndof_reduction_rot = 0
            ntb = int(input_dict['BOUNDCOND'][0])
            ndfmin = int(input_dict['BOUNDCOND'][1])
            if 'COMTRANSROT' in input_dict:
                ncsm = int(input_dict['COMTRANSROT'][0])
                if ncsm < 0 and ntb == 0:
                    # remove rotational COM motion as well - but only if in vacuum!
                    # ncsm < 0 is ignored under PBC
                    if ndfmin != 6:
                        warnings.warn('COMTRANSROT removes both translational and rotational '
                                      'center of mass motion, but BOUNDCOND: NDFMIN = {:d}\n'
                                      'Assuming you know what you are doing and reducing dof by NDFMIN.'
                                      ''.format(ndfmin))
                        system.ndof_reduction_tra = ndfmin
                    else:
                        system.ndof_reduction_tra = 3
                        system.ndof_reduction_rot = 3
                elif ncsm != 0:
                    if ndfmin != 3:
                        warnings.warn('COMTRANSROT removes only translational '
                                      'center of mass motion, but BOUNDCOND: NDFMIN = {:d}\n'
                                      'Assuming you know what you are doing: Reducing dof by NDFMIN.'
                                      ''.format(ndfmin))
                        system.ndof_reduction_tra = ndfmin
                    else:
                        system.ndof_reduction_tra = 3
                elif ndfmin != 0:
                    warnings.warn('COMTRANSROT does not remove '
                                  'center of mass motion, but BOUNDCOND: NDFMIN = {:d}\n'
                                  'Assuming you know what you are doing: Reducing dof by NDFMIN.'
                                  ''.format(ndfmin))
                    system.ndof_reduction_tra = ndfmin
            elif ndfmin != 0:
                warnings.warn('BOUNDCOND: NDFMIN = {:d}, but COMTRANSROT not set.\n'
                              'Assuming you know what you are doing: Reducing dof by NDFMIN.'
                              ''.format(ndfmin))
                system.ndof_reduction_tra = ndfmin

            system.mass = masses

            system.molecule_idx = molecule_idx
            system.nconstraints_per_molecule = nconstraints_per_molec

            system.bonds = all_bonds
            system.constrained_bonds = constrained_bonds

            # thermostat and barostat -> ensemble
            thermostat = False
            temperature = 0
            barostat = False
            pressure = 0
            ambiguous = False

            if 'PRESSURESCALE' in input_dict:
                couple = int(input_dict['PRESSURESCALE'][0])
                scale = int(input_dict['PRESSURESCALE'][1])
                if couple == 2 and scale > 0:
                    barostat = True
                    p = []
                    for n in range(9):
                        p.append(float(input_dict['PRESSURESCALE'][8 + n]))
                    p = np.array(p).reshape((3, 3))
                    p_off = p[~np.eye(p.shape[0], dtype=bool)]
                    p_dia = np.diagonal(p)
                    if not np.allclose(p_off, np.zeros(p_off.shape)):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(non-zero off-diagonal pressure components)')
                        ambiguous = True
                    if not np.all(p_dia == p_dia[0]):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(non-identical diagonal pressure components)')
                        ambiguous = True
                    if scale == 3:  # full anisotropic
                        raise NotImplementedError('TODO: Implement fully anisotropic case.')
                    pressure = p_dia[0]

            if 'MULTIBATH' in input_dict:
                mb_algo = int(input_dict['MULTIBATH'][0])
                if mb_algo == 2:
                    n_nbaths = 2
                else:
                    n_nbaths = 1
                nbaths = int(input_dict['MULTIBATH'][n_nbaths])
                if nbaths == 1:
                    temperature = float(input_dict['MULTIBATH'][n_nbaths + 1])
                    tau = float(input_dict['MULTIBATH'][n_nbaths + 2])
                    if tau > 0:
                        thermostat = True
                elif nbaths > 1:
                    temperatures = np.array([float(t) for t in
                                             input_dict['MULTIBATH'][n_nbaths+1:n_nbaths+nbaths*2:2]])
                    taus = np.array([float(t) for t in
                                     input_dict['MULTIBATH'][n_nbaths+2:n_nbaths+nbaths*2:2]])
                    temperatures = np.array(temperatures[np.where(taus > 0)])
                    if temperatures.size > 0:
                        thermostat = True
                        temperature = temperatures[0]
                        if not np.all(temperatures == temperature):
                            warnings.warn('Ambiguous ensemble definition '
                                          '(baths with different temperature)')
                            ambiguous = True

            if ('STOCHDYN' in input_dict and
               int(input_dict['STOCHDYN'][0]) == 1 and
               int(input_dict['STOCHDYN'][1]) != 0):
                tempsd = float(input_dict['STOCHDYN'][6])
                if not thermostat:
                    thermostat = True
                    temperature = tempsd
                else:
                    if not np.allclose([temperature], [tempsd]):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(SD not at same temperature as MULTIBATH)')
                        ambiguous = True

            volume = None
            if not barostat:
                if 'volume' in observable_dict:
                    volume = np.mean(observable_dict['volume'])
                    if not np.allclose(observable_dict['volume'],
                                       np.ones(observable_dict['volume'].shape) * volume):
                        warnings.warn('No barostat defined, but volume not constant.')
                elif 'volume' in trajectory_dict:
                    volume = np.mean(trajectory_dict['volume'])
                    if not np.allclose(trajectory_dict['volume'],
                                       np.ones(trajectory_dict['volume'].shape) * volume):
                        warnings.warn('No barostat defined, but volume not constant.')

            if not ambiguous:
                if thermostat and barostat:
                    result.ensemble = EnsembleData(
                        ensemble='NPT',
                        natoms=system.natoms,
                        pressure=pressure,
                        temperature=temperature
                    )
                elif thermostat:
                    result.ensemble = EnsembleData(
                        ensemble='NVT',
                        natoms=system.natoms,
                        volume=volume,
                        temperature=temperature
                    )
                else:
                    result.ensemble = EnsembleData(
                        ensemble='NVE',
                        natoms=system.natoms,
                        volume=volume
                    )
                    result.observables['constant_of_motion'] = result.observables['total_energy']
            # end thermostat and barostat -> ensemble
        # end if imd
        # TODO: Have integrator check for existence of constant of motion

        return result
