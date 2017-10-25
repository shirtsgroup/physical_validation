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

from physical_validation.data import parser
from physical_validation.data import simulation_data
from physical_validation.util.gromos_interface import GromosInterface


class GromosParser(parser.Parser):
    """
    GromosParser
    """

    @staticmethod
    def units():
        # Gromos uses kJ/mol
        return simulation_data.UnitData(
            kb=8.314462435405199e-3,
            energy_str='kJ/mol',
            energy_conversion=1.0,
            length_str='nm',
            length_conversion=1.0,
            volume_str='nm^3',
            volume_conversion=1.0,
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
        ensemble: EnsembleData, opitional
            A EnsembleData object
        topology: TopologyData, optional
            A TopologyData object
        tre: str or list, optional
            A string or a list of strings pointing to a .tre file
        cnf: str or list, optional
            A string or a list of strings pointing to a .cnf file 
            (Note: if also trc is given, gro is ignored)
        dt: float, optional
            The time step used in the simulation

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the provided ensemble and
            topology objects as well as the trajectory data found in the
            edr and trr / gro files.

        """
        result = simulation_data.SimulationData()
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
                trv_dict = GromosInterface.read_coords(trc)
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

        result.trajectory = simulation_data.TrajectoryData(
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

            result.observables = simulation_data.ObservableData()
            for key in observables:
                result.observables[key] = observable_dict[key]

        molecules = []
        if top:
            molecules = GromosInterface.read_system_from_top(top)

        input_dict = {}
        natoms = 0
        if imd:
            input_dict = GromosInterface.get_block_fields(imd, unique=True)
            nsolvent_molecules = int(input_dict['SYSTEM'][1])

            # count atoms
            nsolute_atoms = 0
            nsolvent_atoms = 0
            if molecules:
                for molecule in molecules:
                    if molecule['solvent']:
                        molecule['nmolecs'] = nsolvent_molecules
                        nsolvent_atoms += nsolvent_molecules * molecule['natoms']
                    else:
                        nsolute_atoms += molecule['nmolecs'] * molecule['natoms']
            natoms = nsolute_atoms + nsolvent_atoms

            # thermostat and barostat
            thermostat = False
            temperature = 0
            barostat = False
            pressure = 0
            ambiguous = False

            if 'PRESSURESCALE' in input_dict:
                scale = int(input_dict['PRESSURESCALE'][1])
                if scale == 1:
                    barostat = True
                    p = []
                    for n in range(9):
                        p.append(float(input_dict['PRESSURESCALE'][8 + n]))
                    p = np.array(p).reshape((3,3))
                    p_off = p[~np.eye(p.shape[0],dtype=bool)]
                    p_dia = np.diagonal(p)
                    if not np.allclose(p_off, np.zeros(p_off.shape)):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(non-zero off-diagonal pressure components)')
                        ambiguous = True
                    if not np.all(p_dia == p_dia[0]):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(non-identical diagonal pressure components)')
                        ambiguous = True
                    pressure = p_dia[0]

            if 'MULTIBATH' in input_dict:
                nbaths = int(input_dict['MULTIBATH'][1])
                if nbaths == 1:
                    thermostat = True
                    temperature = float(input_dict['MULTIBATH'][2])
                elif nbaths > 1:
                    temperatures = np.array([float(t) for t in input_dict['MULTIBATH'][2:2+nbaths]])
                    if not np.all(temperatures == temperatures[0]):
                        warnings.warn('Ambiguous ensemble definition '
                                      '(baths with different temperature)')
                        ambiguous = True
                    temperature = temperatures[0]
                    thermostat = True

            if ('STOCHDYN' in input_dict and
               int(input_dict['MULTIBATH'][0]) == 1 and
               int(input_dict['MULTIBATH'][1]) != 0):
                tempsd = float(input_dict['MULTIBATH'][6])
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

            ensemble = None
            if not ambiguous:
                if thermostat and barostat:
                    result.ensemble = simulation_data.EnsembleData(
                        ensemble='NPT',
                        natoms=natoms,
                        pressure=pressure,
                        temperature=temperature
                    )
                elif thermostat:
                    result.ensemble = simulation_data.EnsembleData(
                        ensemble='NVT',
                        natoms=natoms,
                        volume=volume,
                        temperature=temperature
                    )
                else:
                    result.ensemble = simulation_data.EnsembleData(
                        ensemble='NVE',
                        natoms=natoms,
                        volume=volume
                    )
                    result.observables['constant_of_motion'] = result.observables['total_energy']

        topology = simulation_data.TopologyData()
        topology.natoms = natoms

        topology.mass =


        if ensemble is not None:
            result.ensemble = ensemble
            if ensemble.ensemble == "NVE":
                self.__gmx_energy_names['constant_of_motion'] = 'Total-Energy'
            else:
                self.__gmx_energy_names['constant_of_motion'] = 'Conserved-En.'

        if topology is not None:
            result.topology = topology

        if edr is not None:
            observable_dict = self.__interface.get_quantities(edr, self.__gmx_energy_names.values())

            # constant volume simulations don't write out the volume in .edr file
            if (observable_dict['Volume'] is None and ensemble is not None and
               ensemble.volume is not None):
                nframes = observable_dict['Pressure'].size
                observable_dict['Volume'] = np.ones(nframes) * ensemble.volume

            result.observables = simulation_data.ObservableData()
            for key, gmxkey in self.__gmx_energy_names.items():
                result.observables[key] = observable_dict[gmxkey]

        if trr is not None:
            if gro is not None:
                warnings.warn('`trr` and `gro` given. Ignoring `gro`.')

            trajectory_dict = self.__interface.read_trr(trr)
            result.trajectory = simulation_data.TrajectoryData(
                trajectory_dict['position'],
                trajectory_dict['velocity'])
        elif gro is not None:
            trajectory_dict = self.__interface.read_gro(gro)
            result.trajectory = simulation_data.TrajectoryData(
                trajectory_dict['position'],
                trajectory_dict['velocity'])

        if dt is not None:
            result.dt = float(dt)

        return result
