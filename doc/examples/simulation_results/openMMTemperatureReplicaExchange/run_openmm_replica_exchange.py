import os

import numpy as np
import openmmtools
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeSampler
from simtk.openmm import XmlSerializer, unit
from simtk.openmm.app.pdbfile import PDBFile

# Run a replica exchange simulation using a stored openmm system.

# Set simulation parameters:
total_simulation_time = 3.0 * unit.nanosecond
simulation_time_step = 5.0 * unit.femtosecond
total_steps = int(np.floor(total_simulation_time / simulation_time_step))
output_directory = "output"
output_data = os.path.join(output_directory, "output.nc")
number_replicas = 6
min_temp = 300.0 * unit.kelvin
max_temp = 500.0 * unit.kelvin
exchange_frequency = 200  # Number of steps between exchange attempts
collision_frequency = 5 / unit.picosecond

# Determine logarithmically spaced temperature list:
T_unit = unit.kelvin
temperature_list = (
    np.logspace(
        np.log10(min_temp.value_in_unit(T_unit)),
        np.log10(max_temp.value_in_unit(T_unit)),
        num=number_replicas,
    )
    * T_unit
)

# Get positions from pdb file:
pdb_path = "initial_structure.pdb"
positions = PDBFile(pdb_path).getPositions()

# Retrieve force field parameters from xml file:
xml_path = "LJ_oligomer_system.xml"
xml_data = open(xml_path).read()
system = XmlSerializer.deserialize(xml_data)

# Determine number of simulation steps to run
simulation_steps = int(np.floor(total_simulation_time / simulation_time_step))

# Determine number of replica exchange attempts
exchange_attempts = int(np.floor(simulation_steps / exchange_frequency))

sampler_states = list()
thermodynamic_states = list()

# Define thermodynamic states
for temperature in temperature_list:
    thermodynamic_state = openmmtools.states.ThermodynamicState(
        system=system, temperature=temperature
    )
    thermodynamic_states.append(thermodynamic_state)
    sampler_states.append(
        openmmtools.states.SamplerState(positions)
    )  # no box vectors, non-periodic system.

# Create and configure simulation object
move = openmmtools.mcmc.LangevinDynamicsMove(
    timestep=simulation_time_step,
    collision_rate=collision_frequency,
    n_steps=exchange_frequency,
    reassign_velocities=False,
)

simulation = ReplicaExchangeSampler(
    mcmc_moves=move,
    number_of_iterations=exchange_attempts,
    replica_mixing_scheme="swap-neighbors",
)

if os.path.exists(output_data):
    os.remove(output_data)

reporter = MultiStateReporter(output_data, checkpoint_interval=1)
simulation.create(thermodynamic_states, sampler_states, reporter)

# Run simulation
simulation.minimize()
simulation.run()
