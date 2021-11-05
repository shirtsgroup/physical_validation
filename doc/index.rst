*****************************
Physical validation reference
*****************************

:code:`physical_validation` is a package aimed at testing results obtained
by molecular dynamics simulations for their physical validity.

.. note:: The physical validation methodology has been described in
   Merz PT, Shirts MR (2018), Testing for physical validity in molecular simulations.
   PLoS ONE 13(9): e0202764. https://doi.org/10.1371/journal.pone.0202764

.. note:: We are always looking to enlarge our set of tests. If you are a
   MD user or developer and have suggestions for physical validity tests
   missing in this package, we would love to hear from you! Please
   consider getting in touch with us via our `github repository`_.

.. toctree::
   userguide
   :maxdepth: 2
   :caption: User guide:

.. toctree::
   examples/kinetic_energy_distribution
   examples/ensemble_check
   examples/kinetic_energy_equipartition
   examples/integrator_validation
   :maxdepth: 2
   :caption: Examples:

.. toctree::
   simulation_data
   :maxdepth: 2
   :caption: Data format and parsers:

.. toctree::
   physical_validation
   physical_validation.data
   physical_validation.util
   :maxdepth: 2
   :caption: Package reference:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`github repository`: https://github.com/shirtsgroup/physical-validation
