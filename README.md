checkensemble
=============

This software allows users to perform statistical test to determine if a given molecular simulation is consistent with the thermodynamic ensemble it is performed in.

Users should cite the JCTC paper: Shirts, "M. R. Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles", J. Chem. Theory Comput., 2013, 9 (2), pp 909–926,  http://dx.doi.org/10.1021/ct300688p


All examples can be run automatically by sourcing 'bash runexamples.sh'

The main distribution consists of the following: 

* checkdist.py
* readmdfiles.py
* analyze-md.py
* analyze-replica.py
* timeseries.py

and python code for toy models for validating the main tools:

* harmonic.py  
* harmonic_pressure.py
* harmonic_mu.py

Run 'analyze-md.py --help' for a listing and explanation of options of
simulation output analysis.
      
Run 'analyze-replica.py --help' for a listing and explanation of
options for temperature replica exchange analysis options

Dependencies: the code requires scipy and numpy.  matplotlib is
required to generate graphs.

For more data, a number of bins closer to 40 may be more appropriate,
but the example data sets are smaller for size reasons, and thus
require smaller bin numbers to get statistically meaningful results.

