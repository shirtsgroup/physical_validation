import numpy as np
from physicalvalidation.data.gromacs_parser import GromacsParser
from physicalvalidation.data.simulation_data import TopologyData, EnsembleData

# SIMULATION DATA CREATION

topo = TopologyData()
topo.natoms = 2700
topo.masses = np.array([15.99940, 1.00800, 1.00800]*900)
topo.nconstraints = 0
topo.ndof_total = 2700*3 - 3
topo.ndof_reduction_tra = 3
topo.ndof_reduction_rot = 0
topo.molecule_idx = np.arange(0, 2700, 3)
topo.nconstraints_per_molecule = np.zeros(900)

NVT_300 = EnsembleData('NVT', natoms=2700, volume=3.01125**3, temperature=300)

parser = GromacsParser(exe='/Users/pascal/Work/gromacs/build/bin/gmx')

nh1_data = parser.get_simulation_data(ensemble=NVT_300, topology=topo,
                                      edr='nh1/water.edr',
                                      gro='nh1/water.gro',
                                      dt=0.0005)

# KINETIC ENERGY VALIDATION

from physicalvalidation import kineticenergy

kineticenergy.check_mb_ensemble(nh1_data,
                                alpha=0.05,
                                verbose=True)

ber1_data = parser.get_simulation_data(ensemble=NVT_300, topology=topo,
                                       edr='ber1/water.edr',
                                       gro='ber1/water.gro',
                                       dt=0.0005)

kineticenergy.check_mb_ensemble(ber1_data,
                                alpha=0.05,
                                verbose=True)

# ENSEMBLE VALIDATION

NVT_310 = EnsembleData('NVT', natoms=2700, volume=3.01125**3, temperature=310)
nh2_data = parser.get_simulation_data(ensemble=NVT_310, topology=topo,
                                      edr='nh2/water.edr',
                                      gro='nh2/water.gro',
                                      dt=0.0005)


from physicalvalidation import ensemble

ensemble.check(nh1_data, nh2_data, total_energy=False)

ber2_data = parser.get_simulation_data(ensemble=NVT_310, topology=topo,
                                       edr='ber2/water.edr',
                                       gro='ber2/water.gro',
                                       dt=0.0005)

ensemble.check(ber1_data, ber2_data, total_energy=False)

# INTEGRATOR CONVERGENCE VALIDATION
nh1_dt_data = parser.get_simulation_data(ensemble=NVT_300, topology=topo,
                                         edr='nh1_dt/water.edr',
                                         gro='nh1_dt/water.gro',
                                         dt=0.00025)

from physicalvalidation import integrator

integrator.convergence([nh1_data, nh1_dt_data], verbose=True, tol=0.1)
