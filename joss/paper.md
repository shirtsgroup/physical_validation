---
title: 'physical_validation: A Python package to assess the physical validity of molecular simulation results'
tags:
  - Python
  - molecular simulation
  - molecular dynamics
  - molecular mechanics
  - statistical mechanics
  - physical validation
authors:
  - name: Pascal T. Merz
    orcid: 0000-0002-7045-8725
    affiliation: 1
  - name: Michael R. Shirts
    orcid: 0000-0003-3249-1097
    affiliation: 1
affiliations:
 - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO 80309, United States of America
   index: 1
date: 10 May 2021
bibliography: paper.bib
---

# Summary

Molecular simulations such as molecular dynamics (MD) and Monte Carlo (MC)
simulations are powerful tools allowing the prediction of experimental
observables in the study of systems such as proteins, membranes, and polymeric
materials.
The quality of predictions based on molecular simulations are depending on the
validity of the underlying physical assumptions.
`physical_validation` allows users of molecular simulation programs to perform
simple tests of physical validity on their systems and setups.
It can also be used by molecular simulation package developers to run
representative test systems during development, increasing code correctness.
The theoretical foundation of the physical validation tests were established
in `[@Merz:2018]`, in which the `physical_validation` package was first
mentioned.

# Statement of need

For a long time, most users of molecular simulation packages were experts that
contributed to the code base themselves or were very familiar with the methodology
used.
Increased popularity of molecular simulation methods have lead to a significantly
increased user base, and to an explosion of available methods.
The simulation packages are faster and more powerful than ever, but still require
expertise to avoid using combination of methods and parameters that could violate
physical assumptions or affect reproducibility.
Unphysical artefacts were reported to significantly affect physical observables
such as the folding of proteins or DNA, the properties of lipid bilayers, the
dynamics of peptides and polymers, or properties of simple liquid (see `[@Merz:2018]`
for further references).
It is a safe assumption that the reported unphysical results are only the tip
of the iceberg, with many more being caught by experienced advisors or reviewers
before publishing, and some being published [by honest mistake, find way to formulate this].

`physical_validation` allows to tackle the problem of robustness in molecular
simulations at two levels.
The first level is the end user.
`physical_validation` allows users to test their simulation results for a number
of deviations from physical assumptions such as the distribution of the kinetic
energy, the equipartition of kinetic energy througout the system, the sampling
of the correct ensemble in the configurational quantities, and the precision of
the integrator.
The combination of these tests allow to cover a wide range of potential
unphysical simulation conditions`[@Merz:2018]`.
This increases the confidence that users can have in their simulation results
independently of code correctness tests provided by the developers of their
simulation package.
The validation tools are explaining their assumptions and conclusions using
comprehensive output and figures, making their use suitable also for users
new to the field of molecular simulations.
Since `physical_validation` is also returning its conslusions in machine-readable
form, it can be included in pipelines allowing results to be tested for
physical validity without user interaction.
The second level are code developers. Unphysical behavior might not only come
from poor or incompatible parameters and models, it might also stem from
coding errors in the simulation programs.
`physical_validation` can be used to regularly run representative simulations
as end-to-end tests in a continuous integration setup, ensuring that code
changes do not introduce bugs that lead to unphysical results.
GROMACS, one of the leading MD packages, has been using `physical_validation`
since 2019 to test every release version with a comprehensive set of
simulations covering all major code paths.

# Acknowledgements

* Grant NIH
* MolSSI
* Can Pervane for helpful discussions in the early stages of the project

Contributions from

* Pascal Merz (279 commits)
* Michael Shirts (6 commits): Original checkensemble, comments on documentation
* Matt Thompson (5 commits): CI
* MatKie (1 commit): Bugfix
* Lenny Fobe (1 commit): Added a tool to CI
* Wei-Tse Hsu (1 commit): Proof-read documentation, still working on additional example
* Nate Abraham (0 commits): Some suggestions for documentation
* Chris Walker: Working on OpenMM replica exchange

# References
