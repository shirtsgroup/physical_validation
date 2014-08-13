#!/usr/bin/python

#=============================================================================================
# example files for reading in MD simulation files and performing
# statistical analyses according to manuscript "Simple tests for
# validity when sampling from thermodynamic ensembles", Michael
# R. Shirts.
# 

# COPYRIGHT NOTICE
#
# Written by Michael R. Shirts <mrshirts@gmail.com>.
#
# Copyright (c) 2012 The University of Virginia. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
# =============================================================================================
#
#===================================================================================================
# IMPORTS
#===================================================================================================

import numpy
import timeseries
import checkensemble 
import optparse, sys
from optparse import OptionParser
import readmdfiles

parser = OptionParser()
parser.add_option("-f", "--files", dest="datafiles",nargs = 2,
                  help="the two files of different temperature for analysis")
parser.add_option("-d", "--directory", dest="datafile_directory",default ='./',
                  help="the directory the data files are is in")
parser.add_option("-k", "--nolikelihood", dest="bMaxLikelihood", action="store_false",default=True,
                  help="Don't run maximum likelihood analysis [default = run this analysis]") 
parser.add_option("-l", "--nolinearfit", dest="bLinearFit", action="store_false",default=True,
                  help="Don't run linear fit analysis [default = run this analysis]") 
parser.add_option("-n", "--nononlinearfit", dest="bNonLinearFit", action="store_false",default=True,
                  help="Don't run linear fit analysis [default = run this analysis]") 
parser.add_option("-t", "--temperature", nargs = 2, dest="T_k", type="float",default=[0,0],
                  help="low and high temperatures, [default = %default]") 
parser.add_option("-p", "--pressure", nargs = 2, dest="P_k", type="float",default=[0,0],
                  help="low and high pressures, [default = %default]") 
parser.add_option("-m", "--chem.potential", nargs = 2, dest="mu_k", type="float",default=[0,0],
                  help="low and high chemical potentials, [default = %default]") 
parser.add_option("-e", "--energytype", dest="type", default="total",
                  help="the type of energy that is being analyzed [default = %default]")
parser.add_option("-b", "--nboot", dest="nboots", type="int",default=200,
                  help="number of bootstrap samples performed [default = %default]")
parser.add_option("-i", "--nbins", dest="nbins",type = "int", default=30,
                  help="number of bins for bootstrapping [default = %default]") 
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",default=False,
                  help="more verbosity")
parser.add_option("-g", "--figurename", dest="figname", default='figure.pdf',
                  help="name for the figure")
parser.add_option("-s", "--seed", dest="seed", type = "int", default=None,
                  help="random seed for bootstrap sampling")
parser.add_option("-c", "--efficiency", nargs = 2, dest="efficiency", type = "float", default=None,
                  help="statistical efficiency to overwrite the calculated statistical efficiency")
parser.add_option("-u", "--useefficiency", dest="useg", type = "string", default='subsample',
                  help= "calculate the efficiency by scaling the uncertainty, or by subsampling the input data")
parser.add_option("--filetype", dest="filetype", type = "string", default='flatfile',
                  help= "specified the type of the file analyzed. options are gromacs .xvg, charmm output, desmond .ene, and flat files")
parser.add_option("--kB", dest="kB", type = "float", default=1.3806488*6.0221415/1000.0,  
                  help="Boltzmann's constant in the applicable units for this sytem")
                  # Boltzmann's constant (kJ/mol/K)

(options, args) = parser.parse_args()

filetypes_supported = ['flatfile','gromacs','charmm','desmond']
if options.datafiles == None:
    print "\nQuitting: No files were input!\n"
    sys.exit()

type = options.type 
onlyE = ['potential', 'kinetic', 'total']
requireV = ['enthalpy', 'volume', 'jointEV'] 
requireN = ['helmholtz', 'number', 'jointEN']
alltypes = onlyE + requireV + requireN

if not (type in alltypes):
    print "type of energy %s isn't defined!" % (type)
    print "Must be one of ", 
    print alltypes
    sys.exit()

if type in onlyE:
    analysis_type = 'dbeta-constV'
elif (type == 'enthalpy'):
    analysis_type = 'dbeta-constP'
elif (type == 'volume'):
    analysis_type = 'dpressure-constB'
elif (type == 'jointEV'):
    analysis_type = 'dbeta-dpressure'
elif (type == 'helmholtz'):
    analysis_type = 'dbeta-constmu'
elif (type == 'number'):
    analysis_type = 'dmu-constB'
elif (type == 'jointEN'):
    analysis_type = 'dbeta-dmu'
else:
    print "analysis type %s not defined: I'll go with total energy" % (type)
    analysis_type = 'dbeta-constV'
if (not(options.useg == 'scale' or options.useg == 'subsample')):
    print "Error: for -u, only options \'scale\' and \'subsample\' allowed"
    sys.exit()

#===================================================================================================
# CONSTANTS
#===================================================================================================

verbose = options.verbose
T_k = numpy.array(options.T_k) #T_k = numpy.array([132.915071475571,137.138128524429])  # temperatures
P_k = numpy.array(options.P_k) #P_k = numpy.array([1.0, 21.0])  # pressures
mu_k = numpy.array(options.mu_k) #mu_k = numpy.array([1.0, 21.0])  # chemical potential
kB = options.kB
names = ['down','up']
nboots = options.nboots
nbins = options.nbins
bMaxLikelihood = options.bMaxLikelihood
bNonLinearFit = options.bNonLinearFit
bLinearFit = options.bLinearFit
figname = options.figname

if not (options.filetype in filetypes_supported):
    print "Error: for -filetype, I currently only know about filetypes",
    print filetypes_supported
    sys.exit()    

if type[0:5] == 'joint':
    bLinearFit = False
    bNonLinearFit = False
    bMaxLikelhood = True
    print "For joint simulations, can only run maximum likelihood, overwriting other options"

if (verbose):
    print "verbosity is %s" % (str(verbose))
    print "Energy type is %s" % (type)
    print "\'%s\' temperature is %f" % (names[0],T_k[0])
    print "\'%s\' temperature is %f" % (names[1],T_k[1])
    if type in requireV: 
       print "\'%s\' pressure is %f" % (names[0],P_k[0])
       print "\'%s\' pressure is %f" % (names[1],P_k[1])
    if type in requireN: 
       print "\'%s\' chemical potential is %f" % (names[0],mu_k[0])
       print "\'%s\' chemical potential is %f" % (names[1],mu_k[1])

    print "Number of bootstraps is %d" % (nboots)
    print "Number of bins (not used for maximum likelihood) is %d" % (nbins)

    if (bMaxLikelihood):
        print "Generating maximum likelihood statistics"
    else:
        print "Not generating maximum likelihood statistics"

    if (bLinearFit):
        print "Generating linear fit statistics"
    else:
        print "Not generating linear fit statistics"

    if (bNonLinearFit):
        print "Generating nonlinear fit statistics"
    else:
        print "Not generating nonlinear fit statistics"


    print "Figures will be named %s" % (figname)    

# Shouldn't need to modify below this for standard usage 
# ------------------------
K = 2;
kJperkcal = 4.184 # unit conversion factor
nm3perA3 = 0.001
N_k = numpy.zeros([K],int) # number of samples at each state

# check just size of all files
N_size = numpy.zeros(K,int) 
filenames = []
for k,T in enumerate(T_k):
    filename = options.datafile_directory + '/' + options.datafiles[k]
    filenames.append(filename)
    print "checking size of \'%s\' file %s..." % (names[k],filenames[k])    
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    N_size[k] = len(lines)

N_max = numpy.max(N_size)
U_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n
V_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # V_kn[k,n] is the volume of the sample k,n
N_kn = numpy.zeros([K,N_max], dtype=int) # U_kn[k,n] is the number of the sample k,n

for k in range(K):
    # Read contents of file into memory.
    print "Reading %s..." % filenames[k]
    infile = open(filenames[k], 'r')
    lines = infile.readlines()
    infile.close()

    if (options.filetype == 'flatfile'): # assumes kJ/mol energies, nm3 volumes
        U_kn[k,:],V_kn[k,:],N_kn[k,:], N_k[k] = readmdfiles.read_flatfile(lines,type,N_max)
        # will need to specify KB for toy models. 
    elif (options.filetype == 'gromacs'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_gromacs(lines,type,N_max)
        N_kn = None
    elif (options.filetype == 'charmm'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_charmm(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
        N_kn = None
    elif (options.filetype == 'desmond'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_desmond(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
        N_kn = None
    else:
        print "The file type %s isn't defined!" % (options.filetype)
        sys.exit()

# compute correlation times for the data
# Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
print "Now determining correlation time"
if (options.efficiency == None):
    g = readmdfiles.getefficiency(N_k,U_kn,V_kn,N_kn,type)
else:
    g = numpy.ones(2);
    print "statistical inefficiencies taken from input options and are %.3f and %.3f steps" % (options.efficiency[0],options.efficiency[1])
if (options.useg == 'subsample'):
    readmdfiles.subsample(N_k,U_kn,V_kn,N_kn,g,type)
else:
    print "statistical efficiencies used to scale the statistical uncertained determined from all data"
figname = options.figname
title = options.figname

checkensemble.ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,P_k=P_k,mu_k=mu_k,U_kn=U_kn,V_kn=V_kn,N_kn=N_kn,title=title,figname=figname,nbins=nbins,reptype='bootstrap',g=g,nboots=nboots,bMaxwell=(type=='kinetic'),bLinearFit=bLinearFit,bNonLinearFit=bNonLinearFit,bMaxLikelihood=bMaxLikelihood,seed=options.seed,kB=kB)
