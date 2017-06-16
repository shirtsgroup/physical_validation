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
gromacs_parser.py
"""
import warnings

import numpy as np

from physical_validation.data import parser
from physical_validation.data import simulation_data
from physical_validation.util.gromacs_interface import GromacsInterface


class GromacsParser(parser.Parser):
    """
    GromacsParser
    """

    @staticmethod
    def units():
        # Gromacs uses kJ/mol
        return simulation_data.UnitData(
            kb=8.314462435405199e-3,
            energy='kJ/mol',
            length='nm',
            volume='nm^3',
            pressure='bar',
            time='ps')

    def __init__(self, exe=None):
        super(GromacsParser, self).__init__()
        self.__interface = GromacsInterface(exe=exe)
        # gmx energy codes
        self.__gmx_energy_names = {'kinetic_energy': 'Kinetic-En.',
                                   'potential_energy': 'Potential',
                                   'total_energy': 'Total-Energy',
                                   'volume': 'Volume',
                                   'pressure': 'Pressure',
                                   'temperature': 'Temperature',
                                   'constant_of_motion': 'Conserved-En.'}

    def get_simulation_data(self,
                            ensemble=None, topology=None,
                            edr=None, trr=None, gro=None,
                            dt=None):
        r"""

        Parameters
        ----------
        ensemble: EnsembleData, opitional
            A EnsembleData object
        topology: TopologyData, optional
            A TopologyData object
        edr: str, optional
            A strint pointing to a .edr file
        trr: str, optional
            A string pointing to a .trr file
        gro: str, optional
            A string pointing to a .gro file (Note: if also trr is given, gro is ignored)
        dt: float, optional
            The time step used in the simulation

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the provided ensemble and
            topology objects as well as the trajectory data found in the
            edr and trr / gro files.

        """
        if list(self.__gmx_energy_names.keys()) != simulation_data.ObservableData.observables():
            warnings.warn('self.__gmx_energy_names.keys() != simulation_data.ObservableData.observables()')

        result = simulation_data.SimulationData()
        result.units = self.units()

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
