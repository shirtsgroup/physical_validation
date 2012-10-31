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
from checkdist import *
import optparse, sys
from optparse import OptionParser
import readmdfiles

parser = OptionParser()
parser.add_option("-f", "--replica_data", dest="metafile", help="prefix of the replica datafiles")
parser.add_option("-k", "--nolikelihood", dest="bMaxLikelihood", action="store_false",default=True,
                  help="Don't run maximum likelihood analysis [default = run this analysis]") 
parser.add_option("-l", "--nolinearfit", dest="bLinearFit", action="store_false",default=True,
                  help="Don't run linear fit analysis [default = run this analysis]") 
parser.add_option("-n", "--nononlinearfit", dest="bNonLinearFit", action="store_false",default=True,
                  help="Don't run linear fit analysis [default = run this analysis]") 
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
parser.add_option("-c", "--efficiency", dest="efficiency", type = "float", default=None,
                  help="statistical efficiency to overwrite the calculated statistical efficiency")
parser.add_option("-u", "--useefficiency", dest="useg", type = "string", default='subsample',
                  help= "calculate the efficiency by scaling the uncertainty, or by subsampling the input data")
parser.add_option("--filetype", dest="filetype", type = "string", default='flatfile',
                  help= "specified the type of the file analyzed. options are gromacs .xvg, charmm output, desmond .ene, and flat files")

(options, args) = parser.parse_args()

if options.metafile == None:
    print "\nQuitting: No files were input!\n"
    sys.exit()

filetypes_supported = ['flatfile','gromacs','charmm','desmond']
onlyE = ['potential', 'kinetic', 'total']
requireV = ['enthalpy', 'volume', 'jointEV'] 
requireN = ['helmholtz', 'number', 'jointEN']
alltypes = onlyE + requireV + requireN

type = options.type

if type in requireN:
    print "Error: replica exchange testing not implemented for grand canonical ensemble yet."

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
# Read metadata.
#===================================================================================================

infile = open(options.metafile, 'r')
lines = infile.readlines()
infile.close()

datafiles = []
T_k = []
P_k = []
K = 0
for line in lines:
    if line[0] != '#':
        elements = line.split()
        K+=1;
        numcol = len(elements)
        printcol = numcol-1
        datafiles.append(elements[0])
        if analysis_type == 'dbeta-constV':
            if (numcol != 2):
                print "Warning! Expecting one temperature entry, getting a different number (%d) of entries!" % (printcol)
            T_k.append(float(elements[1]))
        elif analysis_type == 'dbeta-constP':
            if (numcol != 3):
                print "Warning! Expecting one temperature and on pressure entry, getting a different number (%d) of entries!" % (printcol)
            T_k.append(float(elements[1]))
            P_k.append(float(elements[2]))
        elif analysis_type == 'dpressure-constB':
            if (numcol != 2):
                print "Warning! Expecting one pressure entry, getting a different number (%d) of entries!" % (printcol)
            P_k.append(float(elements[1]))
        elif analysis_type == 'dbeta-dpressure':
            if (numcol != 3):
                print "Warning! Expecting one temperature and on pressure entry, getting a different number (%d) of entries!" % (printcol)
            T_k.append(float(elements[1]))
            P_k.append(float(elements[2]))

T_k = numpy.array(T_k)
P_k = numpy.array(P_k)

#===================================================================================================
# CONSTANTS
#===================================================================================================

verbose = options.verbose
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
    for k in range(K):
        print "%dth temperature is %f" % (k,T_k[k])
        if ((type == 'volume') or (type == 'enthalpy') or (type == 'jointPV')):                
            print "%dth pressure is %f" % (k,P_k[k])
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
kB = 1.3806488*6.0221415/1000.0  # Boltzmann's constant (kJ/mol/K)   - gromacs default
kJperkcal = 4.184 # unit conversion factor
nm3perA3 = 0.001
N_k = numpy.zeros([K],int) # number of samples at each state

# check just size of all files
N_size = numpy.zeros(K,int) 
filenames = []
for k in range(K):
    filename = datafiles[k]
    filenames.append(filename)
    print "checking size of file \#%d named %s..." % (k+1,filenames[k])    
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    N_size[k] = len(lines)

N_max = numpy.max(N_size)
U_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n
V_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n
N_kn = None  # replica exchange doesn't support different chemical potentials yet.

for k in range(K):
    # Read contents of file into memory.
    print "Reading %s..." % filenames[k]
    infile = open(filenames[k], 'r')
    lines = infile.readlines()
    infile.close()

    if (options.filetype == 'flatfile'): # assumes kJ/mol energies, nm3 volumes
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_flatfile(lines,type,N_max)
    elif (options.filetype == 'gromacs'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_gromacs(lines,type,N_max)
    elif (options.filetype == 'charmm'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_charmm(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
    elif (options.filetype == 'desmond'):
        U_kn[k,:],V_kn[k,:],N_k[k] = readmdfiles.read_desmond(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
    else:
        print "The file type %s isn't defined!" % (options.filetype)
        sys.exit()

# compute correlation times for the data
# Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
print "Now determining correlation time"
if (options.efficiency == None):
    g = readmdfiles.getefficiency(N_k,U_kn,V_kn,N_kn,type)
else:
    g = options.efficiency*numpy.ones(K)
    print "statistical inefficiency taken from input options is %f" % (options.efficiency)
if (options.useg == 'subsample'):
    readmdfiles.subsample(N_k,U_kn,V_kn,N_kn,g,type)
else:
    print "statistical efficiencies used to scale the statistical uncertained determined from all data"
figname = options.figname
title = options.figname

for k in range(K-1):
    twoN = numpy.array([N_k[k],N_k[k+1]])
    if (type in onlyE) or (type == 'enthalpy') or (type == 'jointEV'):
        twoT = numpy.array([T_k[k],T_k[k+1]])
    else:
        twoT = None
    if type in requireV:
        twoP = numpy.array([P_k[k],P_k[k+1]])
    else:
        twoP = None
    twoU = U_kn[k:k+2,:]
    twoV = V_kn[k:k+2,:]
    ProbabilityAnalysis(twoN,type=analysis_type,T_k=twoT,P_k=twoP,U_kn=twoU,V_kn=twoV,nbins=nbins,
                        reptype='bootstrap',g=g,nboots=nboots,bMaxwell=(type=='kinetic'),figname='replica',bLinearFit=bLinearFit,bNonLinearFit=bNonLinearFit,bMaxLikelihood=bMaxLikelihood,seed=options.seed)

    # now, we construct a graph with all of the lines. We could write
    # the probability analysis to do it, but better to do new specially designed plot here. 

    
