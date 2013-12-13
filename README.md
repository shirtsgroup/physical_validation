#checkensemble
#=============
#
#This software allows users to perform statistical test to determine if a given molecular simulation is consistent with the thermodynamic ensemble it is performed in.

# This README file can be sourced as 'bash README.txt' to automatically run all examples.
# The main distribution consists of the following: 
#
# checkdist.py
# readmdfiles.py
# analyze-md.py
# analyze-replica.py
# timeseries.py
#
# and python code for toy models for validating the main tools.
#
# harmonic.py  
# harmonic_pressure.py
# harmonic_mu.py
#
# Run 'analyze-md.py --help' for a listing and explanation of options of simulation output analysis.
#      
# Run 'analyze-replica.py --help' for a listing and explanation of options for 
#     temperature replica exchange analysis options
#
# Dependencies: the code requires scipy and numpy.  matplotlib is required to generate graphs.
#
# For more data, a number of bins closer to 40 may be more appropriate, but the example 
# data sets are smaller for size reasons, and thus require smaller bin numbers to get significant ):
#
##################
#Analysis of gromacs files for vrescale thermostat:

python analyze-md.py --filetype gromacs -d example_data -f vrescale.argon.down.xvg vrescale.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 20 -v -g vrescale -e total -s 12345 > vrescale.txt

#Produces figures vrescale_linear.pdf and vrescale_nonlinear.pdf

##################

#Analysis of gromacs files for berendsen thermostat:

python analyze-md.py --filetype gromacs -d example_data -f berendsen.argon.down.xvg berendsen.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 20 -v -g berendsen -e total -s 12345 > berendsen.txt

#Produces figures berendsen_linear.pdf and berendsen_nonlinear.pdf

##################

#Analysis of gromacs files for Parrinell-Rahman barostat - tests enthaply distribution

python analyze-md.py --filetype gromacs -d example_data -f prH.argon.down.xvg prH.argon.up.xvg -t 121.43116059 128.56883941 -p 90 90 -b 200 -i 20 -v -g enthalpy -e enthalpy -s 12345 > enthalpy.txt

#Produces figures enthalpy_linear.pdf and enthalpy_nonlinear.pdf

##################

#Analysis of gromacs files for Parrinell-Rahman barostat - tests volume distribution

python analyze-md.py --filetype gromacs -d example_data -f prV.argon.down.xvg prV.argon.up.xvg -t 125 125 -p 30 150 -b 200 -i 20 -v -g volume -e volume -s 12345 > volume.txt
#Produces figures volume_linear.pdf and volume_nonlinear.pdf

##################

#Analysis of gromacs files for Parrinell-Rahman barostat - test joint energy/volume distribution

python analyze-md.py --filetype gromacs -d example_data -f prEV.argon.down.xvg prEV.argon.up.xvg -t 121.43116059 128.56883941  -p 30 150 -b 200 -i 20 -v -g energyvolume -e jointEV -s 12345 > energyvolume.txt
#Produces figures energyvolume_linear.pdf and energyvolume_nonlinear.pdf

##################

#Analysis of charmm files of water run with the Nose-Hoover thermostat - demonstrating alternate file types.
#No bootstraps are run because there are not enough samples to avoid NaN's in all of the 
#histograms.  charmm output files have relatively few samples per MB!  If only maximum
#likelihood is run, bootstrap sampling can be performed.

python analyze-md.py --filetype charmm -d example_data -f charmm_300.out charmm_303.out -t 300 303 -b 0 -i 20 -g charmm -e total > charmm.txt
 
##################

#Analysis of desmond files of run with the Nose-Hoover thermostat.  Currently dummy files -- run at the same temperature,
#and too little data, but demonstrates the file type support.

python analyze-md.py --filetype desmond -d example_data -f Desmond_test_1.ene Desmond_test_2.ene -t 298 298 -b 0 -c 1 1 -i 3 -g desmond -e total > desmond.txt

##################

#Analysis of flat files that other simulation data can be written in; just repurposing of earlier gromacs files. (vrescale.argon,prV.argon,prEV.argon)

python analyze-md.py --filetype flatfile -d example_data -f flat.argon.E.down.txt flat.argon.E.up.txt -t 132.915071475571 137.138128524429 -b 0 -i 20 -g flat_E -e total -v > flat.E.output.txt
python analyze-md.py --filetype flatfile -d example_data -f flat.argon.V.down.txt flat.argon.V.up.txt -t 125 125 -p 30 150 -b 0 -i 20 -g flat_V -e volume -v > flat.V.output.txt
python analyze-md.py --filetype flatfile -d example_data -f flat.argon.EV.down.txt flat.argon.EV.up.txt -t 121.43116059 128.56883941 -p 30 150 -b 0 -i 20 -e jointEV -v > flat.EV.output.txt

##################

#Analysis of harmonic oscillators:

python harmonic.py -t 0.5 1.5 -g harmonic_energy -s 12345 > harmonic_energy.txt

#Produces figures harmonic_energy_linear.pdf and harmonic_energy_nonlinear.pdf

##################

#harmonic oscillators with some noise, creating deviations from the fit.

python harmonic.py -t 0.5 1.5 -o 0.1 -g harmonic_energy_noise > harmonic_energy_noise.txt

#Produces figures harmonic_energy_noise_linear.pdf and harmonic_energy_noise_nonlinear.pdf
##################

#Analysis of harmonic oscillators with work functions:

python harmonic_pressure.py -t 0.5 1.5 -p 1.0 1.0 -e enthalpy -g harmonic_enthalpy -s 12345 > harmonic_enthalpy.txt
#Produces figures harmonic_enthalpy_linear.pdf and harmonic_enthalpy_nonlinear.pdf

python harmonic_pressure.py -t 1.0 1.0 -p 0.7 1.3 -e volume -g harmonic_volume -s 12345 > harmonic_volume.txt
#Produces figures harmonic_volume_linear.pdf and harmonic_volume_nonlinear.pdf

python harmonic_pressure.py -t 1.7 1.3 -p 0.8 1.2 -e jointEV -s 12345 > harmonic_jointEV.txt

##################

#Analysis of changing numbers of harmonic oscillators (Grand Canonical Ensemble)

python harmonic_mu.py -t 0.5 0.5 -m -0.5 -0.4 -e number -g harmonic_number -s 12345 > harmonic_number.txt
#Produces figures harmonic_number_linear.pdf and harmonic_number_nonlinear.pdf

python harmonic_mu.py -t 0.5 0.6 -m -0.5 -0.5 -e helmholtz -g harmonic_helmholtz -s 12345 > harmonic_helmholtz.txt
#Produces figures harmonic_helmholtz_linear.pdf and harmonic_helmholtz_nonlinear.pdf

python harmonic_mu.py -t 0.4 0.5 -m -0.5 -0.4 -e jointEN -g harmonic_jointEN -s 12345 > harmonic_jointEN

##################

#Testing flat file validation for Grand Canonical Monte Carlo simulations

python analyze-md.py --filetype flatfile -d example_data -f flat.GC.N.down.txt flat.GC.N.up.txt -t 0.5 0.5 -m -0.5 -0.4 --kB 1 -e number -g flat_N -s 12345 > flatfile.N.output.txt
#Produces figures flat_N_linear.pdf and flat_N_nonlinear.pdf

python analyze-md.py --filetype flatfile -d example_data -f flat.GC.A.down.txt flat.GC.A.up.txt -t 0.5 0.6 -m -0.5 -0.5 --kB 1 -e helmholtz -g flat_A -s 12345 > flatfile.A.output.txt
#Produces figures flat_A_linear.pdf and flat_A_nonlinear.pdf

python analyze-md.py --filetype flatfile -d example_data -f flat.GC.EN.down.txt flat.GC.EN.up.txt analyze_md.py -t 0.4 0.5 -m -0.5 -0.4 --kB 1 -e jointEN -s 12345 > flatfile.EN.output.txt

##################

#Replica exchange analysis.
python analyze-replica.py --filetype gromacs -f example_data/replica_files.txt -b 40 -l -n -c 5 -i 20 -v -g vrescale_kinetic -e kinetic -s 12345 > replica_total_energy.txt

python analyze-replica.py --filetype gromacs -f example_data/replica_files.txt -b 40 -l -n -c 5 -i 20 -v -g vrescale_total -e total -s 12345 > replica_kinetic_energy.txt

