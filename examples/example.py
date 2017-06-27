from physical_validation.data.gromacs_parser import GromacsParser

# Replace this with path to exectuable if gmx is not in your path
parser = GromacsParser(exe='gmx')

nh1_data = parser.get_simulation_data(mdp='nh1/water.mdp',
                                      top='nh1/water.top',
                                      edr='nh1/water.edr',
                                      gro='nh1/water.gro')

# KINETIC ENERGY VALIDATION

from physical_validation import kinetic_energy

kinetic_energy.mb_ensemble(nh1_data,
                           alpha=0.05,
                           verbose=True)

ber1_data = parser.get_simulation_data(mdp='ber1/water.mdp',
                                       top='ber1/water.top',
                                       edr='ber1/water.edr',
                                       gro='ber1/water.gro')

kinetic_energy.mb_ensemble(ber1_data,
                           alpha=0.05,
                           verbose=True)

# ENSEMBLE VALIDATION

nh2_data = parser.get_simulation_data(mdp='nh2/water.mdp',
                                      top='nh2/water.top',
                                      edr='nh2/water.edr',
                                      gro='nh2/water.gro')


from physical_validation import ensemble

ensemble.check(nh1_data, nh2_data, total_energy=False)

ber2_data = parser.get_simulation_data(mdp='ber2/water.mdp',
                                       top='ber2/water.top',
                                       edr='ber2/water.edr',
                                       gro='ber2/water.gro')

ensemble.check(ber1_data, ber2_data, total_energy=False)

# INTEGRATOR CONVERGENCE VALIDATION
nh1_dt_data = parser.get_simulation_data(mdp='nh1_dt/water.mdp',
                                         top='nh1_dt/water.top',
                                         edr='nh1_dt/water.edr',
                                         gro='nh1_dt/water.gro')

from physical_validation import integrator

integrator.convergence([nh1_data, nh1_dt_data], verbose=True, tol=0.1)
