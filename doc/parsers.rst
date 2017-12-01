.. _doc_parsers:

:class:`.SimulationData` objects and parsers
=============================================

The data of simulations to be validated are best represented by objects
of the  :class:`.SimulationData` type. While lower-level functions accepting
bare arrays and numbers are available, the  :class:`.SimulationData` objects
combine ease-of-use and higher stability in form of input testing.

The  :class:`.SimulationData` objects are consisting of information about the
simulation and the system. This information is collected in objects of different
classes, namely

* :obj:`.SimulationData.units` of type :class:`.UnitData`:
  Information on the units used by the simulation program.
* :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`:
  Information on the sampled ensemble.
* :obj:`.SimulationData.system` of type :class:`.SystemData`:
  Information on the system (atoms, molecules, constraints, etc.).
* :obj:`.SimulationData.observables` of type :class:`.ObservableData`:
  Trajectories of observables along the simulation.
* :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`:
  Position / velocity / force trajectories along the simulation.
* :obj:`.SimulationData.dt` of type `float`:
  The time step at which the simulation was performed.

The :class:`.SimulationData` objects can either be constructed
directly from arrays and numbers, or (partially) automatically via parsers.

Package-specific parsers
------------------------

Package-specific parsers return :class:`.SimulationData` objects via the
:func:`.Parser.get_simulation_data` function by reading the output files
of the corresponding MD program.

GROMACS: :class:`.GromacsParser`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :class:`.GromacsParser` takes the GROMACS input files `mdp` (run options)
and `top` (topology file) to read the details about the system, the ensemble
and the time step. The observable trajectory is extracted from an `edr`
(binary energy trajectory), while the position and velocity trajectory can
be read either from a `trr` (binary trajectory) or a `gro` (ascii trajectory)
file. The constructor optionally takes the path to a gromacs binary as well
as the path to the topology library as inputs. The first is necessary to
extract information from binary files (using `gmx energy` and `gmx dump`),
while becomes necessary if the `top` file contains `#include` statements
which usually rely on GROMACS environment variables.

Example usage:
::

   import physical_validation as pv

   parser = pv.data.GromacsParser(exe='~/bin/gromacs/bin/gmx',
                                  includepath='~/bin/gromacs/share/gromacs/top')

   res = parser.get_simulation_data(
        mdp='mdout.mdp',
        top='system.top',
        gro='system.gro',
        edr='system.edr'
   )

Flatfile parser
---------------

For MD packages not supported by the package-specific parsers, there is the
possibility to create the :class:`.SimulationData` objects via the
:class:`.FlatfileParser`. This parser fills the
:obj:`.SimulationData.trajectory` object via 3-dimensional ascii files
containing the position and velocity trajectories, and the
:obj:`.SimulationData.observables` via 1-dimensional ascii files containing
the trajectories for the observables of interest. As the details on the
units, the simulated system and the sampled ensemble can not easily be read
from such files, this information has to be provided by the user by passing
objects of the respective data structures. See
:func:`.FlatfileParser.get_simulation_data` for more details on the
:class:`.SimulationData` creation via the flat file parser, and
:ref:`simulationdata_details` for details on which test requires which
information.


Example usage, system of 900 water molecules in GROMACS units simulated in
NVT (note that this example leaves some fields in :class:`.SystemData`
empty, as well as the trajectory of some observables and the position and
velocities):
::

   import physical_validation as pv

   parser = pv.data.FlatfileParser()

   system = pv.data.SystemData(
       natoms=900*3,
       nconstraints=900*3,
       ndof_reduction_tra=3,
       ndof_reduction_rot=0
   )

   units = pv.data.UnitData(
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
       time_conversion=1.0
   )

   ensemble = pv.data.EnsembleData(
       ensemble='NVT',
       natoms=900*3,
       volume=3.01125**3,
       temperature=298.15
   )

   res = parser.get_simulation_data(
       units=units, ensemble=ensemble, system=system,
       kinetic_ene_file='kinetic.dat',
       potential_ene_file='potential.dat',
       total_ene_file='total.dat'
   )


Create :class:`.SimulationData` objects from python data
--------------------------------------------------------

As an alternative to the different parsers described above,
:class:`.SimulationData` objects can of course also be created by calling
the constructor directly. The :class:`.SimulationData` constructor thereby
simply takes the objects listed at the beginning of this section. All these
objects can also be set later during the lifetime of the object - a
:class:`.SimulationData` object can hence be initialized empty and filled
at a later point. The objects contained in :class:`.SimulationData` objects
are explained in details in :ref:`simulationdata_details`.

Example usage, system of 900 water molecules in GROMACS units simulated in
NVT (note that this example leaves some fields in :class:`.SystemData`
empty, as well as the trajectory of some observables and the position and
velocities):
::

   import physical_validation as pv

   system = pv.data.SystemData(
       natoms=900*3,
       nconstraints=900*3,
       ndof_reduction_tra=3,
       ndof_reduction_rot=0
   )

   units = pv.data.UnitData(
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
       time_conversion=1.0
   )

   ensemble = pv.data.EnsembleData(
       ensemble='NVT',
       natoms=900*3,
       volume=3.01125**3,
       temperature=298.15
   )

   # This snippet is assuming that kin_ene, pot_ene and tot_ene are lists
   # or numpy arrays filled with the kinetic, potential and total energy
   # of a simulation run. These might be obtained, e.g., from the python
   # API of a simulation code, or from other python-based analysis tools.

   observables = pv.data.ObservableData(
       kinetic_energy = kin_ene,
       potential_energy = pot_ene,
       total_energy = tot_ene
   )

   res = pv.data.SimulationData(
       units=units, ensemble=ensemble,
       system=system, observables=observables
   )


.. _simulationdata_details:

Data contained in :class:`.SimulationData` objects
==================================================

Units: :obj:`.SimulationData.units` of type :class:`.UnitData`
--------------------------------------------------------------
Attributes:

* :attr:`.UnitData.kb`, `float`
* :attr:`.UnitData.energy_conversion`, `float`
* :attr:`.UnitData.length_conversion`, `float`
* :attr:`.UnitData.volume_conversion`, `float`
* :attr:`.UnitData.temperature_conversion`, `float`
* :attr:`.UnitData.pressure_conversion`, `float`
* :attr:`.UnitData.time_conversion`, `float`
* :attr:`.UnitData.energy_str`, `str`
* :attr:`.UnitData.length_str`, `str`
* :attr:`.UnitData.volume_str`, `str`
* :attr:`.UnitData.temperature_str`, `str`
* :attr:`.UnitData.pressure_str`, `str`
* :attr:`.UnitData.time_str`, `str`

The information about units consists of different parts:
* The value of kB in the used energy units,
* the conversion factor to GROMACS units (kJ/mol, nm, nm^3, K, bar, ps), and
* the name of the units (energy_str, length_str, volume_str, temperature_str, pressure_str, time_str).
The names are only used for output (console printing and plotting), and are optional.
The conversion factors and kB are, on the other hand, used in computations and need
to be given.

Needed by

  * :func:`physical_validation.ensemble.check`
  * :func:`physical_validation.ensemble.estimate_interval`
  * :func:`physical_validation.kinetic_energy.mb_ensemble`, only

    - :attr:`.UnitData.kb`

Ensemble: :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`
------------------------------------------------------------------------
Attributes:

* :attr:`.EnsembleData.ensemble`, `str`
* :attr:`.EnsembleData.natoms`, `int`
* :attr:`.EnsembleData.mu`, `float`
* :attr:`.EnsembleData.volume`, `float`
* :attr:`.EnsembleData.pressure`, `float`
* :attr:`.EnsembleData.energy`, `float`
* :attr:`.EnsembleData.temperature`, `float`

The ensemble is a string indicating the thermodynamical ensemble a simulation was
performed in, and is any of 'NVE', 'NVT', 'NPT', 'muVT'.

Depending on the ensemble, EnsembleData then holds additional information defining
the ensemble, such as the number of particles N, the chemical potential mu, the
volume V, the pressure P, the constant energy E or the temperature T. While any
of these additional information are optional, most of them are needed by certain
tests, such that not fully defining the ensemble results in warnings. The notable
exception to this rule is the constant energy E for NVE, which is not needed
by any test and can hence be omitted without raising a warning.

Needed by
  * :func:`physical_validation.kinetic_energy.mb_ensemble`
  * :func:`physical_validation.ensemble.check`

System: :obj:`.SimulationData.system` of type :class:`.SystemData`
------------------------------------------------------------------
Attributes:

    * :attr:`.SimulationData.natoms`, the total number of atoms in the system;
      e.g. for a system containing 100 water molecules: `.SimulationData.natoms = 300`
    * :attr:`.SimulationData.nconstraints`, the total number of constraints in the
      system, not including the global translational and rotational constraints
      (see next two attributes); e.g. for a system containing 100 *rigid* water molecules:
      `.SimulationData.nconstraints = 300`
    * :attr:`.SimulationData.ndof_reduction_tra`, global reduction of translational
      degrees of freedom (e.g. due to constraining the center of mass of the system)
    * :attr:`.SimulationData.ndof_reduction_rot`, global reduction of rotational
      degrees of freedom (e.g. due to constraining the center of mass of the system)
    * :attr:`.SimulationData.mass`, a list of the mass of every atom in the system;
      e.g. for a single water molecule: `.SimulationData.mass = [15.9994, 1.008, 1.008]`
    * :attr:`.SimulationData.molecule_idx`, a list with the indices first atoms of every
      molecule (this assumes that the atoms are sorted by molecule); e.g. for a system
      containing 3 water molecules: `.SimulationData.molecule_idx = [0, 3, 6]`
    * :attr:`.SimulationData.nconstraints_per_molecule`, a list with the number of
      constraints in every molecule; e.g. for a system containing 3 *rigid* water
      molecules: `.SimulationData.nconstraints_per_molecule = [3, 3, 3]`
    * :attr:`.SimulationData.bonds`, a list containing all bonds in the system;
      e.g. for a system containing 3 water molecules:
      `.SimulationData.bonds = [[0, 1], [0, 2], [3, 4], [3, 5], [6, 7], [6, 8]]`
    * :attr:`.SimulationData.constrained_bonds`, a list containing only the constrained
      bonds in the system, must be a subset of `.SimulationData.bonds` (and equal, if
      all bonds are constrained).

.. todo:: Currently, there is some redundancy in the attributes listed above. The
   :attr:`.SimulationData.bonds` and :attr:`.SimulationData.constrained_bonds` are
   reserved for future use - included already in the information about the system,
   but not yet used by any tests included in the currently published package. In a
   future version, the :class:`.SystemData` should be streamlined to make the object
   initialization easier.

Needed by

  * :func:`physical_validation.kinetic_energy.mb_ensemble`, partially:

    - :attr:`.SystemData.natoms`,
    - :attr:`.SystemData.nconstraints`,
    - :attr:`.SystemData.ndof_reduction_tra`,
    - :attr:`.SystemData.ndof_reduction_rot`

  * :func:`physical_validation.kinetic_energy.equipartition`, all attributes except
    :attr:`.SimulationData.bonds` and :attr:`.SimulationData.constrained_bonds`.

Observables: :obj:`.SimulationData.observables` of type :class:`.ObservableData`
--------------------------------------------------------------------------------
Attributes:

  * :attr:`.ObservableData.kinetic_energy`, the kinetic energy trajectory (nframes x 1),
    also accessible via `.ObservableData['kinetic_energy']`
  * :attr:`.ObservableData.potential_energy`, the potential energy trajectory (nframes x 1),
    also accessible via `.ObservableData['potential_energy']`
  * :attr:`.ObservableData.total_energy`, the total energy trajectory (nframes x 1),
    also accessible via `.ObservableData['total_energy']`
  * :attr:`.ObservableData.volume`, the volume trajectory (nframes x 1),
    also accessible via `.ObservableData['volume']`
  * :attr:`.ObservableData.pressure` the pressure trajectory (nframes x 1),
    also accessible via `.ObservableData['pressure']`
  * :attr:`.ObservableData.temperature` the temperature trajectory (nframes x 1),
    also accessible via `.ObservableData['temperature']`
  * :attr:`.ObservableData.constant_of_motion` the constant of motion trajectory (nframes x 1),
    also accessible via `.ObservableData['constant_of_motion']`

Needed by

  * :func:`physical_validation.kinetic_energy.mb_ensemble`

    - :attr:`.ObservableData.kinetic_energy`

  * :func:`physical_validation.ensemble.check`

    - :attr:`.ObservableData.total_energy`, or
    - :attr:`.ObservableData.potential_energy`,
    - :attr:`.ObservableData.volume` (for NPT)

  * :func:`physical_validation.integrator.convergence`

    - :attr:`.ObservableData.constant_of_motion`

Atom trajectories: :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`
-------------------------------------------------------------------------------------
Attributes:

  * :attr:`.TrajectoryData.position`, the position trajectory (nframes x natoms x 3),
    also accessible via `.TrajectoryData['position']`
  * :attr:`.TrajectoryData.velocity`, the velocity trajectory (nframes x natoms x 3),
    also accessible via `.TrajectoryData['velocity']`

Needed by

  * :func:`physical_validation.kinetic_energy.equipartition`


Time step: :obj:`.SimulationData.dt` of type `float`
----------------------------------------------------
The timestep used during the simulation run, a single `float` value.

Needed by

  * :func:`physical_validation.integrator.convergence`
