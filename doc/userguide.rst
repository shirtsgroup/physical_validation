Introduction
============

Advances in recent years have made molecular dynamics (MD) simulations a
powerful tool in molecular-level research, allowing the prediction of
experimental observables in the study of systems such as proteins, drug
targets or membranes. The quality of any prediction based on MD results
will, however, strongly depend on the validity of underlying physical
assumptions.

This package is intended to help detect (sometimes hard-to-spot)
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
parsers for the outputs of several popular MD programs are provided

Second, it can be integrated in MD code testing environments. While
unphysical behavior can be due to poor or incompatible choices of
parameters by the user, it can also originate in coding errors
within the program. Physical validation tests can be integrated in the
code-checking mechanism of MD software packages to facilitate the
detection of such bugs. The :code:`physical_validation` package is currently
used in the automated code-testing facility of the GROMACS software
package, ensuring that every major releases passes a number of physical
sanity checks performed on selected representative systems before
shipping.

.. note:: The physical validation tests have been described in [Merz2018]_.

.. note:: We are always looking to enlarge our set of tests. If you are a
   MD user or developer and have suggestions for physical validity tests
   missing in this package, we would love to hear from you! Please
   consider getting in touch with us via our `github repository`_.


Installation
============

pip
---
The most recent release of `physical_validation` can be installed from `PyPI`_
via :code:`pip`
::

   pip install physical_validation

conda
-----
The most recent release of `physical_validation` can also be installed using
:code:`conda`
::

   conda install -c conda-forge physical_validation

Development version
-------------------

The latest version is available on our `github repository`_. You can install
it via :code:`pip`
::

   pip install git+https://github.com/shirtsgroup/physical_validation.git


Simulation data
===============

The data of simulations to be validated are represented by objects
of the :class:`.SimulationData` type.

The :class:`.SimulationData` objects contain information about the simulation
and the system. This information is collected in objects of different
classes, namely:

* :obj:`.SimulationData.units` of type :class:`.UnitData`:
  Information on the units used by the simulation program.
* :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`:
  Information on the sampled ensemble. This includes the temperature, pressure, and chemical potential,
  with specific requirements depending on the ensemble specified.
* :obj:`.SimulationData.system` of type :class:`.SystemData`:
  Information on the system (number of atoms, molecules, constraints, etc.).
* :obj:`.SimulationData.observables` of type :class:`.ObservableData`:
  Trajectories of observables along the simulation, such as energy or volume. 
* :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`:
  Position / velocity / force trajectories along the simulation.
* :obj:`.SimulationData.dt` of type :code:`float`:
  The time step at which the simulation was performed.

The different physical validation tests do not all require all data to be
able to run. Each :code:`physical_validation` function checks whether the required
information was provided, and raises an error if the information is
insufficient. :ref:`simulationdata_details` lists by which tests the single
members of :class:`.SimulationData` are required.

The :class:`.SimulationData` objects can either be constructed directly
from arrays and numbers, or (partially) automatically via parsers.
The preferred way to populate :class:`.SimulationData` objects is by
assigning its sub-objects explicitly with data obtained from the simulation
package. Many simulation packages have a well-defined Python interface which
allows to read observable, position and velocity trajectories into Python data
structures. The remaining information, such as details on the simulated
ensemble or the molecular system, is usually rather easy to fill in by hand.
The examples in this documentation follow this model.

Please see :ref:`doc_simulation_data` for more details on the
:class:`.SimulationData` type and on how to create objects from results
obtained from different simulation packages.


Kinetic energy validation
=========================
Kinetic energy validation includes testing the likelihood of a trajectory
to originate from the theoretically expected gamma distribution and
validating the temperature equipartition between groups of degrees
of freedom. For details on the employed algorithms, please check the
respective function documentations.

For both the full distribution test and the equipartition test, a strict
and a non-strict version are available. They are triggered using the
:code:`strict=[True|False]` keyword. The strict version does a full distribution
similarity analysis using the Kolmogorov-Smirnov (K-S) test. The K-S test
returns a p-value indicating the likelihood that the sample originates from
the expected distribution. Its sensitivity
increases with increasing sample size, and can flag even the smallest deviations
from the expected distribution at large sample sizes. When developing or
implementing new temperature control algorithms in a controlled testing
environment which keeps errors from other sources negligible, such a high
sensibility is desirable. In other
applications, however, a deviation insignificant in comparison with
other sources of inaccuracies might be enough to flag long simulation
trajectories of large systems as not having a gamma distribution. For
example, deviations from the desired kinetic energy distribution that
are smaller in magnitude than other well-controlled approximations, such as
the interaction cutoff or the treatment of bond constraints, might be enough
to flag large samples as not being properly distributed.

As an alternative to the strict test, the :code:`physical_validation` suite offers
the non-strict version. In this case, the mean and the standard deviation of
the sample are calculated and compared to the expected values. To make the
test easily interpretable, two distinct temperatures :math:`T_\mu` and
:math:`T_\sigma` are estimated from the kinetic energy distribution. They represent the
temperature at which the sample mean and standard would be physically expected.
An error estimate computed via bootstrapping of the provided kinetic energy samples is given for each of the
temperatures, giving information on the statistical significance of the results.

For more details about the difference between the strict test and non-strict test, please
see :func:`physical_validation.kinetic_energy.distribution`.

Full system distribution validation
-----------------------------------
Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.kinetic_energy.distribution`

Example
~~~~~~~
`Kinetic energy distribution example`_

.. _`Kinetic energy distribution example`: examples/kinetic_energy_distribution.ipynb

Equipartition validation
------------------------
Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.kinetic_energy.equipartition`

Example
~~~~~~~
`Kinetic energy equipartition example`_

.. _`Kinetic energy equipartition example`: examples/kinetic_energy_equipartition.ipynb


Ensemble validation
===================
As the distribution of configurational quantities like the potential
energy :math:`U`, the volume :math:`V` or (for the grand and semigrand canonical ensembles) 
the number of each species :math:`N_i` are in general not known analytically, testing the likelihood
of a trajectory sampling a given ensemble is less straightforward than
for the kinetic energy. However, generally, the *ratio* of the probability
distribution between samplings of the same system generated at different state
points (e.g. simulations run at at different temperatures or different pressures) is exactly known for each ensemble
[Merz2018]_, [Shirts2013]_.
Providing two simulations at different state points therefore allows a
validation of the sampled ensemble.

Note that the ensemble validation function is automatically inferring the
correct test based on the simulation input data (such as temperature and pressure) that are given as input.

Choice of the state points
--------------------------
As the above ensemble tests require two simulations at distinct
state points, the choice of interval between the two points is an
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

When setting :code:`verbosity >= 1` in :func:`physical_validation.ensemble.check`, the
routine is printing an estimate for the optimal spacing based on the distributions
provided. Additionally, :func:`physical_validation.ensemble.estimate_interval`
calculates the estimate given a single simulation result. This can be used to determine
at which state point a simulation should be repeated in order to efficiently check
its sampled ensemble.

Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.ensemble.check`

:func:`physical_validation.ensemble.estimate_interval`

Example
~~~~~~~
`Ensemble validation example`_

.. _`Ensemble validation example`: examples/ensemble_check.ipynb


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

Example
-------
`Integrator convergence example`_

.. _`Integrator convergence example`: examples/integrator_validation.ipynb


.. _`PyPI`: https://pypi.org/project/physical_validation

.. _`github repository`: https://github.com/shirtsgroup/physical_validation

.. [Merz2018] Merz PT, Shirts MR (2018)
   "Testing for physical validity in molecular simulations",
   PLOS ONE 13(9): e0202764.
   https://doi.org/10.1371/journal.pone.0202764

.. [Shirts2013] Shirts, M.R.
   "Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles",
   J. Chem. Theory Comput., 2013, 9 (2), pp 909–926,
   http://dx.doi.org/10.1021/ct300688p
