r"""
gromacs_parser.py
"""
import warnings

import numpy as np

from physicalvalidation.data import parser
from physicalvalidation.data import simulation_data
from physicalvalidation.util.gromacs_interface import GromacsInterface


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

    def __init__(self, exe=None, dp=None):
        super(GromacsParser, self).__init__()
        self.__interface = GromacsInterface(exe=exe, dp=dp)
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
                            edr=None, trr=None, gro=None):

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

        return result
