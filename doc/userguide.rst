Introduction
============

Advances in recent years have made molecular dynamics (MD) simulations a
powerful tool in molecular-level research, allowing the prediction of
experimental observables in the study of systems such as proteins, drug
targets or membranes. The quality of any prediction based on MD results
targets or membranes. The quality of any prediction based on MD results
will, however, strongly depend on the validity of underlying physical
assumptions.

This package is intended to help detecting (sometimes hard-to-spot)
unphysical behavior of simulations, which may have statistically important
influence on their results. It is part of a two-fold approach to
increase the robustness of molecular simulations.

First, it empowers users of MD programs to test the physical validity on
their respective systems and setups. The tests range from simple
post-processing analysis to more involved tests requiring additional
simulations. These tests can significantly increase the
reliability of MD simulations by catching a number of common simulation
errors violating physical assumptions, such as non-conservative
integrators, deviations from the specified Boltzmann ensemble, or lack of ergodicity
between degrees of freedom. To make usage as easy as possible,
parsers for several popular MD program output formats are provided.

Second, it can be integrated in MD code testing environments. While
unphysical behavior can be due to poor or incompatible choices of
parameters by the user, it can also originate in coding errors
within the program. Physical validation tests can be integrated in the
code-checking mechanism of MD software packages to facilitate the
detection of such bugs. The `physical-validation` package is currently
used in the automated code-testing facility of the GROMACS software
package, ensuring that every major releases passes a number of physical
sanity checks performed on selected representative systems before
shipping.

.. note:: We are always looking to enlarge our set of tests. If you are a
   MD user or developer and have suggestions for physical validity tests
   missing in this package, we would love to hear from you! Please
   consider getting in touch with us via our `github page`_.


Installation
============

The latest version is available at the `github page`_. Once you downloaded
or cloned the repository, simply type
::

   python3 setup.py install

while being in the root directory of the downloaded package.


`examples/` folder
==================

The folder `examples/` contains examples (simulation results and analysis
scripts) which are used in the following to introduce the functionality of
the package. More specifically, the subfolders of `examples/` contain

* Folder `water_ensemble/`: GROMACS result files of water systems simulated with
  different thermostating and barostating algorithms. Specifically, the
  subfolder `be/` contains NVT results obtained using a Berendsen thermostat,
  the subfolder `vr/` contains NVT results obtained using the GROMACS v-rescale
  thermostat, while the subfolders `be_pr/` and `vr_pr/` contain the corresponding
  NPT results obtained by adding a Parinello-Rahman barostat to the systems. Each
  system was simulated at two different state points (different temperatures, and
  possibly pressures), stored in subfolders `base/` and `ensemble_1/` under the
  respective system folders. Each of these folders then contain the following
  files:

  - `start.gro`: the starting configuration, containing 900 three-site water molecules
  - `system.top`: the topology of the system of water molecules
  - `system.mdp`: the GROMACS input file
  - `mdout.mdp`: the GROMACS input file obtained from `gmx grompp`
  - `system.gro`: the end configuration
  - `system.edr`: the resulting (binary) energy file

* Folder `argon_integrator/`: GROMACS result files of argon systems simulated with
  different van-der-Waals cutoff schemes. `none/` contains the results of
  simulations performed with vdW interactions cut off at 1 nm without any
  correction of the discontinuity at the cut-off distance. `shift/` contains
  the results of simulations performed with vdW interactions cut off at 1 nm
  and a potential shifted by a constant value such that it reaches zero at
  the cut-off distance. `switch/` contains the results of simulations performed
  with vdW interactions cut off at 1 nm where not only the potential, but also
  the forces were altered in a way to have both the forces and the potential
  smoothly reaching zero at the cut-off distance. For each choice of interaction
  function, the same simulation was repeated five times, with each new simulation
  halving the integration time step compared to the previous one.


Simulation data
===============

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

To facilitate the use of the physical validation suite, simulation results
generated by selected software packages can be automatically created by
:class:`.Parser` objects. The currently supported MD packages are:

* GROMACS: :class:`.GromacsParser`

Package-specific parsers are subclasses of :class:`.Parser`, and need to
redefine the :func:`.Parser.get_simulation_data` returning a
:class:`.SimulationData` object.

For generic input, the flat file parser :class:`.FlatfileParser` allows to
create a :class:`.SimulationData` object from files containing trajectories
of observables and the positions and velocities of the atoms in the system.
It requires, however, to fill information on the units, the ensemble and the
system by hand. Furthermore, it is of course possible to fill all attributes
of :class:`.SimulationData` by hand, e.g. when starting from data stored
in python arrays rather than in files.

Please see :ref:`doc_parsers` for more details on the :class:`.SimulationData`
type and the available parsers.

.. note:: We are looking to enlarge the collection of parsers to make the
   use of the package as convenient as possible for as many users as
   possible. If your MD program of choice is not supported (yet), please
   consider either writing your own parser and contribute it by creating
   a pull request on the project's `github page`_, or contacting us to
   let us know about your needs, and we can coordinate about introducing
   the appropriate parser.

.. _example_sec_1:

Examples
--------
To illustrate the creation of SimulationData, we will look at the first part
of the analysis script `ana_water.py` located in the `examples/water_ensemble/`
folder. First, after some necessary import and definitions, the GROMACS
parser is created:
::

   import physical_validation as pv
   import os

   systems = ['vr', 'be', 'vr_pr', 'be_pr']

   # change this to fit your GROMACS installation
   parser = pv.data.GromacsParser(exe='~/bin/gromacs/bin/gmx',
                                  includepath='~/bin/gromacs/share/gromacs/top')

Having the parser readily available, actually reading in the simulation data is a
one-line command that is easily included in a loop for the different systems
of interest:
::

   for sys in systems:
       print('### Analyzing system ' + sys)
       print('## Reading lower temperature result')
       dir_low = os.path.join(sys, 'base')
       res_low = parser.get_simulation_data(
           mdp=os.path.join(dir_low, 'mdout.mdp'),
           top=os.path.join(dir_low, 'system.top'),
           gro=os.path.join(dir_low, 'system.gro'),
           edr=os.path.join(dir_low, 'system.edr')
       )
       print('## Reading high temperature result')
       dir_high = os.path.join(sys, 'ensemble_1')
       res_high = parser.get_simulation_data(
           mdp=os.path.join(dir_high, 'mdout.mdp'),
           top=os.path.join(dir_high, 'system.top'),
           gro=os.path.join(dir_high, 'system.gro'),
           edr=os.path.join(dir_high, 'system.edr')
       )


.. _ke_val_full:

Kinetic energy validation (full system)
=======================================
Kinetic energy validation includes testing the likelihood of a trajectory
to originate from the theoretically expected gamma distribution and
validating the temperature equipartition between groups of degrees
of freedom. For details on the employed algorithms, please check the
respective function documentations. For an example of the equipartition
test, see :ref:`ke_val_equi`.

For both the full distribution test and the equipartition test, a strict
and a non-strict version are available. They are trigged using the
`strict={True|False}` keyword. The strict version does a full distribution
similarity analysis using the Kolmogorov-Smirnov (K-S) test. The K-S test
returns a p-value indicating the likelihood that the sample originates from
the expected distribution. Its sensitivity
increases with increasing sample size, and can flag even smallest deviations
from the expected distribution at large sample sizes. When developing or
implementing new temperature control algorithms in a controlled test
environment keeping errors from other sources negligible, a very high
sensibility is certainly desirable. In many other, real-world
applications, however, a deviation insignificant in comparison with
other sources of inaccuracies might be enough to flag long simulation
trajectories of large systems as not gamma distributed. For
example, deviations from the desired kinetic energy distribution which
are smaller in magnitude than other well-controlled approximations such as
the interaction cutoff or the treatment of bond constraints might be enough
to flag large samples as not being properly distributed.

As an alternative to the strict test, the `physical_validation` suite offers
the non-strict version. In this case, the mean and the standard deviation of
the sample are calculated and compared to the expected values. To make the
test easily interpretable, a distinct temperature :math:`T_\mu` and
:math:`T_\sigma` is associated to each of the two moments. They represent the
temperature at which the sample mean and width would be physically expected.
An error estimate computed via bootstrapping is given for each of the
temperatures, giving information on the statistical significance of the results.

Functions
---------
*Distribution validation:*
:func:`physical_validation.kinetic_energy.distribution`

Examples
--------
With the data structures created in :ref:`example_sec_1` (`res_low` and
`res_high`), the kinetic energy ensemble of each simulated state point
can be validated via the strict test as follows:
::

   print('\n## Validating kinetic energy distribution (strict)')
   print('# Low T:')
   pv.kinetic_energy.distribution(res_low, verbosity=2, strict=True,
                                  filename=sysplot + '_low_mb')
   print('# High T:')
   pv.kinetic_energy.distribution(res_high, verbosity=2, strict=True,
                                  filename=sysplot + '_high_mb')

This will plot the sampled distribution along with its analytical counterpart,
and print out the result of the analysis. For the NVT simulation
using the v-rescale algorithm (folder `vr/base`), as an example, the result will
indicate that the null-hypothesis that the energy is Maxwell-Boltzmann distributed
cannot be rejected:
::

   After equilibration, decorrelation and tail pruning, 96.02% (4802 frames) of original Kinetic energy remain.
   Kinetic energy distribution check (strict)
   Kolmogorov-Smirnov test result: p = 0.901109
   Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed

On the other hand, the NVT simulation using the Berendsen algorithm will show
a dramatically different picture:
::

   After equilibration, decorrelation and tail pruning, 95.96% (4799 frames) of original Kinetic energy remain.
   Kinetic energy distribution check (strict)
   Kolmogorov-Smirnov test result: p = 2.00541e-17
   Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed

The non-strict test confirms this finding, and actually gives a hint to the
reason (that can easily be confirmed by looking at the plotted distribution):
The distribution sampled by the Berendsen algorithm is significantly too narrow.
::

   After equilibration, decorrelation and tail pruning, 84.64% (4233 frames) of original Kinetic energy remain.
   Kinetic energy distribution check (non-strict)
   Analytical distribution (T=298.15 K):
    * mu: 6689.47 kJ/mol
    * sigma: 128.77 kJ/mol
   Trajectory:
    * mu: 6692.54 +- 1.99 kJ/mol
      T(mu) = 298.29 +- 0.09 K
    * sigma: 128.38 +- 1.28 kJ/mol
      T(sigma) = 297.23 +- 2.95 K

::

   After equilibration, decorrelation and tail pruning, 95.96% (4799 frames) of original Kinetic energy remain.
   Kinetic energy distribution check (non-strict)
   Analytical distribution (T=298.15 K):
    * mu: 6689.47 kJ/mol
    * sigma: 128.77 kJ/mol
   Trajectory:
    * mu: 6690.21 +- 1.44 kJ/mol
      T(mu) = 298.18 +- 0.06 K
    * sigma: 98.81 +- 1.03 kJ/mol
      T(sigma) = 228.78 +- 2.37 K


Ensemble validation
===================
As the distribution of configurational quantities like the potential
energy :math:`U`, the volume :math:`V` or (for the grand and semigrand canonical ensembles) 
the number of each species are in general not known analytically, testing the likelihood
of a trajectory sampling a given ensemble is less straightforward than
for the kinetic energy. However, generally, the ratio of the probability
distribution between samplings of the same ensemble at different state
points (e.g. at different temperatures, different pressures) is known
[Shirts2013]_.
Providing two simulations at different state points therefore allows a
validation of the sampled ensemble.

Note that the ensemble validation function is automatically inferring the
correct test based on the simulation that are given as input.

.. [Shirts2013] Shirts, M.R.
   "Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles",
   J. Chem. Theory Comput., 2013, 9 (2), pp 909–926,
   http://dx.doi.org/10.1021/ct300688p

Functions
---------
:func:`physical_validation.ensemble.check`

Examples
--------
Still using the data structures created in :ref:`example_sec_1` (`res_low` and
`res_high`), the generated ensemble of the potential energy can now be validated,
to check whether a similar trend as for the kinetic energy can be observed. The
relevant line of code reads
::

   print('\n## Validating ensemble')
   quantiles = pv.ensemble.check(res_low, res_high, quiet=False,
                                 screen=False, filename=sysplot + '_ensemble')

The ensemble validation function uses the two simulation results at lower and
higher state point to calculate the ratio of the energy distributions and
compare this ratio to the analytical expectation. Under NVT conditions, the
validation will always compare simulations at different temperatures.
Under NPT conditions, there are three possible tests, each requiring
different simulations. If the two simulations were performed at
different temperatures, then the distribution of the instantaneous
enthalpy :math:`U + PV` is tested.  If the two simulations were
performed at different pressures, then the distribution of :math:`V`
is tested. If simulations were performed at both different
temperatures and pressures, then test of the joint distribution of
:math:`U` and :math:`V` is performed (as in our example below).

The relevant result
from these calculations is the deviation from the analytical expectation,
reported in terms of the number of standard deviations (quantiles) the result
is off. The maximum-likelihood analysis of the NPT simulations performed with
the v-rescale thermostat and the Parrinello-Rahman barostat reads
::

   ==================================================
   Maximum Likelihood Analysis (analytical error)
   ==================================================
   Free energy
       521.66705 +/- 13.00986
   Estimated slope                  |  True slope
       0.012963  +/- 0.000301       |  0.013091
       (0.42 quantiles from true slope)
       -2.232868 +/- 0.157176       |  -2.349681
       (0.74 quantiles from true slope)
   Estimated dT                     |  True dT
       9.9    +/- 0.2               |  10.0
   Estimated dP                     |  True dP
       93.5   +/- 6.6               |  98.3
   ==================================================
   Calculated slope is (0.4, 0.7) quantiles from the true slope

This corresponds to a near-perfect agreement with the analytical expectation,
suggesting that the ensemble sampled by the potential energy and the volume
is very close to a the desired NPT ensemble.

Performing the same analysis with the NPT simulations using the Berendsen
thermostat and the Parrinello-Rahman barostat leads to a significantly
different result:
::

   ==================================================
   Maximum Likelihood Analysis (analytical error)
   ==================================================
   Free energy
       900.85325 +/- 25.81833
   Estimated slope                  |  True slope
       0.023477  +/- 0.000653       |  0.013091
       (15.90 quantiles from true slope)
       -2.446500 +/- 0.226693       |  -2.349681
       (0.43 quantiles from true slope)
   Estimated dT                     |  True dT
       17.9   +/- 0.5               |  10.0
   Estimated dP                     |  True dP
       102.4  +/- 9.5               |  98.3
   ==================================================
   Calculated slope is (15.9, 0.4) quantiles from the true slope

This result indicates that using Berendsen thermostat does not only not
generate the proper distribution of the kinetic energy, but does also
affect the ratio of potential energy distribution at different
temperatures. The pressure distribution, governed by the
Parrinello-Rahman barostat, on the other hand, does not seem affected.

Note that for both the NVT and the NPT ensemble, the test involving
different temperatures can also be performed using the total energy
:math:`U + K` (NVT) or :math:`U + PV + K` (NPT). This option can be
enabled using the `total_energy = True` flag of the
:func:`physical_validation.ensemble.check` function, which is disabled
by default. As the kinetic energy can be checked separately (see above),
using the total energy will in general not give any additional insights
and might mask errors in the other energy terms.

Support for grand and semigrand canonical ensembles, validating the
distribution of $N$ and $U$ or composition will be provided soon; in
the meantime, this functionality can still be found in the
checkensemble_ repository.

Choice of the state points
--------------------------
As the ensemble tests presented above require two simulations at distinct
state points, the choice of interval between the two points becomes an
important question. Choosing two state points too far apart will result
in poor or zero overlap between the distributions, leading to very noisy
results (due to sample errors in the tails) or a breakdown of the method,
respectively. Choosing two state points very close to each others, on the
other hand, makes it difficult to distinguish the slope from statistical
error in the samples.

A rule of thumb states [Shirts2013]_ that the maximal efficiency of the
method is reached when the distance between the peaks of the distributions
are roughly equal to the sum of their standard deviations. For most systems
with the exception of extremely small or very cold systems, it is reasonable
to assume that the difference in standard deviations between the state points
will be negligable. This leads to two ways of calculating the intervals:

*Using calculated standard deviations*: Given a simulation at one state point,
the standard deviation of the distributions can be calculated numerically. The
suggested intervals are then given by

* :math:`\Delta T = 2 k_B T^2 / \sigma_E`, where :math:`\sigma_E` is the standard
  deviation of the energy distribution used in the test (potential energy, enthalpy,
  or total energy).
* :math:`\Delta P = 2 k_B T / \sigma_V`, where :math:`\sigma_V` is the standard
  deviation of the volume distribution.

*Using physical observables*: The standard deviations can also be estimated using
physical observables such as the heat capacity and the compressibility. The
suggested intervals are then given by:

* :math:`\Delta T = T (2 k_B / C_V)^{1/2}` (NVT), or
  :math:`\Delta T = T (2 k_B / C_P)^{1/2}` (NPT), where :math:`C_V` and :math:`C_P`
  denote the isochoric and the isobaric heat capacities, respectively.
* :math:`\Delta P = (2 k_B T / V \kappa_T)`, where :math:`\kappa_T` denotes the
  isothermal compressibility.

When setting `verbosity >= 1` in :func:`physical_validation.ensemble.check`, the
routine is printing an estimate for the optimal spacing based on the distributions
provided. Additionaly, :func:`physical_validation.ensemble.estimate_interval`
calculates the estimate given a single simulation result. This can be used to determine
at which state point a simulation should be repeated in order to efficiently check
its sampled ensemble.


.. _ke_val_equi:

Kinetic energy validation (equipartition)
=========================================

Functions
---------
*Equipartition validation:*
:func:`physical_validation.kinetic_energy.equipartition`

.. todo:: Equipartition example

Integrator Validation
=====================
A symplectic integrator can be shown to conserve a constant of motion
(such as the energy in a microcanonical simulation) up to a fluctuation
that is quadratic in time step chosen. Comparing two or more
constant-of-motion trajectories realized using different time steps (but
otherwise unchanged simulation parameters) allows a check of the
symplecticity of the integration. Note that lack of symplecticity does not
necessarily imply an error in the integration algorithm, it can also hint
at physical violations in other parts of the model, such as non-continuous
potential functions, imprecise handling of constraints, etc.

Functions
---------
:func:`physical_validation.integrator.convergence`

Examples
--------
To demonstrate the integration validation, we will use the results in
folder `argon_integrator`, and the corresponding analysis script
`ana_argon.py` located in that folder. As described above, this folder
contains the result of an argon system simulated with different cut-off
schemes of the van-der-Waals interactions, two of which include a
discontinuity of the forces (in subfolder `shift/`) or even
both the forces and the potential (in subfolder `none/`). These small
discontinuities are not unlike bugs that could be present in an interaction
calculation, and will therefore be used to demonstrate the use of the
integrator convergence validation to detect errors in MD codes.

The first lines of the script `ana_argon.py` are very similar to the
previously discussed script `ana_water.py`, setting up the necessary
prerequisites and reading the results using the parser. The actual test
is then called as
::

   pv.integrator.convergence(res, verbose=True,
                             filename=sysplot)

where `res` is a list of :class:`.SimulationData` objects of identical
simulations performed at different integrator time steps.

The final output of the script `ana_argon.py` reads
::

   ### Analyzing system none
   ## Reading results
   ## Validating integrator convergence
   -----------------------------------------------------------------
           dt        avg       rmsd      slope         ratio
                                                     dt^2       rmsd
   -----------------------------------------------------------------
        0.004   -4749.12   3.66e-01   2.86e-04         --         --
        0.002   -4749.27   3.72e-01   2.77e-04       4.00       0.99
        0.001   -4749.26   3.34e-01   3.47e-04       4.00       1.11
       0.0005   -4749.23   3.37e-01   3.33e-04       4.00       0.99
      0.00025   -4749.23   3.45e-01   2.54e-04       4.00       0.98
   -----------------------------------------------------------------

   ### Analyzing system shift
   ## Reading results
   ## Validating integrator convergence
   -----------------------------------------------------------------
           dt        avg       rmsd      slope         ratio
                                                     dt^2       rmsd
   -----------------------------------------------------------------
        0.004   -4491.08   1.63e-02  -1.76e-07         --         --
        0.002   -4491.24   4.51e-03  -1.98e-06       4.00       3.62
        0.001   -4491.24   1.36e-03  -2.55e-06       4.00       3.31
       0.0005   -4491.21   2.83e-04  -2.46e-07       4.00       4.81
      0.00025   -4491.19   1.20e-04   2.96e-07       4.00       2.35
   -----------------------------------------------------------------

   ### Analyzing system switch
   ## Reading results
   ## Validating integrator convergence
   -----------------------------------------------------------------
           dt        avg       rmsd      slope         ratio
                                                     dt^2       rmsd
   -----------------------------------------------------------------
        0.004   -4335.09   1.69e-02   5.54e-07         --         --
        0.002   -4335.25   4.37e-03  -4.87e-07       4.00       3.87
        0.001   -4335.24   1.09e-03  -3.81e-08       4.00       4.02
       0.0005   -4335.22   2.77e-04  -2.67e-08       4.00       3.93
      0.00025   -4335.20   6.90e-05  -9.41e-09       4.00       4.02
   -----------------------------------------------------------------

The outputs of the function are the time step, the average value of the
constant of motion, and its RMSD during the simulation. The fourth
column gives the measured slope of the constant of motion - a large
value here would indicate a strong drift and hence a problem in the
integrator. Even without strong drift, as in the current situation, a
large deviation in the ratio between the rmsd values compared to the
ratio between the time step will indicates some error in the integrator.
The reason for a failure of this test might not always be intuitively clear,
as many components play into the integrator convergence - the integrator
algorithm itself, but also the interaction function (e.g. non-continuous
cut-off) or the numerical precision of the floating point operations.

In the examples presented here, the integrator convergence validation
shows a high sensibility towards the incontinuities describes above. In
the case with discontinuous potential and forces, the constant of motion
shows practically no dependence on the time step. But also with the
shifted (and hence continuous) potential, the large fluctuations around
the expected convergence indicate a problem in the calculation. Ensuring
continuity in the forces allows, on the other hand, to massively reduce
these fluctuations.

.. _`github page`: https://github.com/shirtsgroup/physical-validation

.. _checkensemble: https://github.com/shirtsgroup/checkensemble
