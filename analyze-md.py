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

def read_flatfile(lines,type,N_max):

    # assumes kJ/mol energies, nm^3 volumes 
    # allocate space

    U_n = numpy.zeros([N_max], dtype=numpy.float64) # U_n[k,n] is the energy of the sample n
    V_n = numpy.zeros([N_max], dtype=numpy.float64) # V_n[k,n] is the volume of the sample n
    N = 0

    # we assume energy is first, then volume
    for line in lines:
        if (line[0] != '#'):   # in flat file format, anything that starts with a hash is an ignored comment
        
            elements = line.split()
            numcol = len(elements)
            if (numcol == 0):
                print "Error: No data for data point %d" % (N) 
                sys.exit()
            if (type == 'enthalpy') or (type == 'jointEV'):
                if (numcol != 2):
                    print "Error: asking for enthalpy or jointEV test but %d data points provided instead of 2" % (numcol)
                    sys.exit()
            else:
                if (numcol != 1):
                    print "Error: asking for energy or volume test but %d data point provided instead of 1" % (numcol)
                    sys.exit()

            if (type != 'volume'):
                U_n[N] = float(elements[0])
            if (type == 'volume'): 
                V_n[N] = float(elements[0])
            if (type == 'enthalpy') or (type == 'jointEV'):
                V_n[N] = float(elements[1])
            N += 1 

    return U_n,V_n,N        

def read_gromacs(lines,type,N_max):

    # allocate space
    U_n = numpy.zeros([N_max], dtype=numpy.float64) # U_n[k,n] is the energy of the sample n
    V_n = numpy.zeros([N_max], dtype=numpy.float64) # V_n[k,n] is the volume of the sample n

    N = 0
    ematch = False
    vmatch = False
    for line in lines:
        # Split line into elements.
        if (line[0:3] == '@ s'):
            elements = line.split()
            whichcol = int((elements[1])[1:])+1   # figure out which column it is
            if (type == 'potential'):
                if (elements[3] == "\"Potential\""):
                    ecol = whichcol
                    ematch = True
            if (type == 'total') or (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                if (elements[3] == "\"Total"):
                    comp = elements[3] + ' ' + elements[4]
                    if (comp == "\"Total Energy\""):
                        ecol = whichcol
                        ematch = True
            if (type == 'kinetic'):
                if (elements[3] == "\"Kinetic"):
                    comp = elements[3] + ' ' + elements[4]
                    if (comp == "\"Kinetic En.\""):
                        ecol = whichcol
                        ematch = True
            if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                if (elements[3] == "\"Volume\""):
                    vcol = whichcol
                    vmatch = True

        if ((line[0] != '#') and (line[0] != '@')):
                
           elements = line.split()
           # what is the time of the sample
           if (type != 'volume'):
               energy = float(elements[ecol])
               U_n[N] = energy   
           if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
               volume = float(elements[vcol])
               V_n[N] = volume   
           N += 1 

    return U_n,V_n,N      

def read_charmm(lines,type,N_max):

    # allocate space
    U_n = numpy.zeros([N_max], dtype=numpy.float64) # U_n[k,n] is the energy of the sample n
    V_n = numpy.zeros([N_max], dtype=numpy.float64) # V_n[k,n] is the volume of the sample n

    N = 0
    ematch = False
    vmatch = False
    for line in lines:
        elements = line.split()
        if (line[0:4] == 'DYNA'):
            if (line[0:8] == 'DYNA DYN'):
                if (line[8:9] == ':'):
                    for i,e in enumerate(elements):
                        if (type == 'kinetic'):
                            if (e[0:4] == 'TOTK'):
                                ecol = i-1
                                ematch = True
                        if (type == 'potential'):
                            if (e[0:4] == 'ENER'):
                                ecol = i-1
                                ematch = True
                        if (type == 'total') or (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                            if (e[0:4] == 'TOTE'):
                                ecol = i-1
                                ematch = True
            elif (line[0:5] == 'DYNA>'):
                U_n[N] = float(elements[ecol])   
                if (type != 'volume'):
                    N += 1   # we count here unless volume is the only variable

            if (line[0:10] == 'DYNA PRESS'):
                if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                    if (line[10:11] == ':'):
                        for i,e in enumerate(elements):
                            if (e[0:4] == 'VOLU'):
                                vcol = i
                                vmatch = True
                    elif (line[10:11] == '>'):                
                        V_n[N] = float(elements[vcol])   
                        if (type == 'volume'):
                            N += 1   # we only count here when volume is the only variable
    return U_n,V_n,N      

def read_desmond(lines,type,N_max):

    # reads desmond .ene files

    # allocate space
    U_n = numpy.zeros([N_max], dtype=numpy.float64) # U_n[k,n] is the energy of the sample n
    V_n = numpy.zeros([N_max], dtype=numpy.float64) # V_n[k,n] is the volume of the sample n

    N = 0
    ematch = False
    vmatch = False
    for line in lines:
        # Split line into elements.
        if (line[0:13] == '# 0:time (ps)'):    # this line tells us what column is what
            elements = line.split()
            for e in (elements):
                if (e[0] != '(') and e[0] != '#':
                    (num,val) = e.split(':')
                    if (type == 'potential'):
                        if val == 'E_p':
                            ecol = int(num)
                            ematch = True
                    if (type == 'total') or (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                        if val == 'E':
                            ecol = int(num)
                            ematch = True
                    if (type == 'kinetic'):
                        if val == 'E_k':
                            ecol = int(num)
                            ematch = True
                    if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
                        if val == 'V':
                            vcol = int(num)
                            vmatch = True

        if ((line[0] != '#') and (line != '\n')):
                
           elements = line.split()
           # what is the time of the sample
           if (type != 'volume'):
               energy = float(elements[ecol])
               U_n[N] = energy   
           if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
               volume = float(elements[vcol])
               V_n[N] = volume   
           N += 1 

    return U_n,V_n,N      

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
parser.add_option("-t", "--temperature", nargs = 2, dest="T_k", type="float",default=[295,305],
                  help="low and high temperatures, [default = %default]") 
parser.add_option("-p", "--pressure", nargs = 2, dest="P_k", type="float",default=[1,21],
                  help="low and high pressures, [default = %default]") 
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
parser.add_option("--filetype", dest="filetype", type = "string", default='gromacs',
                  help= "specified the type of the file analyzed. options are gromacs .xvg, charmm output, desmond .ene, and flat files")

(options, args) = parser.parse_args()

type = options.type
if (type != 'kinetic') and (type != 'potential') and (type != 'total') and (type != 'enthalpy') and (type != 'volume') and (type != 'jointEV'):
    print "type of energy %s isn't defined!" % (type)
    sys.exit()

if ((type == 'kinetic') or (type == 'potential') or (type == 'total')):
    analysis_type = 'dbeta-constV'
elif (type == 'enthalpy'):
    analysis_type = 'dbeta-constP'
elif (type == 'volume'):
    analysis_type = 'dpressure-constB'
elif (type == 'jointEV'):
    analysis_type = 'dbeta-dpressure'
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
names = ['down','up']
type = options.type # 'potential', 'kinetic', 'total', 'enthalpy', 'volume', 'jointEV'
nboots = options.nboots
nbins = options.nbins
bMaxLikelihood = options.bMaxLikelihood
bNonLinearFit = options.bNonLinearFit
bLinearFit = options.bLinearFit
figname = options.figname

if (not(options.filetype == 'gromacs' or options.filetype == 'charmm' or options.filetype == 'desmond' or options.filetype == 'flatfile')):
    print "Error: for -filetype, I currently only know about filetypes \'flatfile\' \'gromacs\', \'charmm\', and \'desmond\'."
    sys.exit()    

if (type == 'jointEV'):
    bLinearFit = False
    bNonLinearFit = False
    bMaxLikelhood = True
    print "For type \'JointPV\' can only run maximum likelihood, overwriting other options"

if (verbose):
    print "verbosity is %s" % (str(verbose))
    print "Energy type is %s" % (type)
    print "\'%s\' temperature is %f" % (names[0],T_k[0])
    print "\'%s\' temperature is %f" % (names[1],T_k[1])
    if ((type == 'volume') or (type == 'enthalpy') or (type == 'jointPV')):                
       print "\'%s\' pressure is %f" % (names[0],P_k[0])
       print "\'%s\' pressure is %f" % (names[1],P_k[1])
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
kB = 1.3806488*6.0221415/1000.0  # Boltzmann's constant (kJ/mol/K)   - gromacs default
kJperkcal = 4.184 # unit conversion factor
nm3perA3 = 0.001
N_k = numpy.zeros([K],int) # number of samples at each state

# check just size of all files
N_size = numpy.zeros(K,int) 
filenames = []
for k,T in enumerate(T_k):
    filename = options.datafile_directory + '/' + options.datafiles[k]
    filenames.append(filename)
    print "checking size of \'%s\' temperature file %s..." % (names[k],filenames[k])    
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    N_size[k] = len(lines)

N_max = numpy.max(N_size)
U_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n
V_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n

for k,T in enumerate(T_k):

    # Read contents of file into memory.
    print "Reading %s..." % filenames[k]
    infile = open(filenames[k], 'r')
    lines = infile.readlines()
    infile.close()

    if (options.filetype == 'flatfile'): # assumes kJ/mol energies, nm3 volumes
        U_kn[k,:],V_kn[k,:],N_k[k] = read_flatfile(lines,type,N_max)
    elif (options.filetype == 'gromacs'):
        U_kn[k,:],V_kn[k,:],N_k[k] = read_gromacs(lines,type,N_max)
    elif (options.filetype == 'charmm'):
        U_kn[k,:],V_kn[k,:],N_k[k] = read_charmm(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
    elif (options.filetype == 'desmond'):
        U_kn[k,:],V_kn[k,:],N_k[k] = read_desmond(lines,type,N_max)
        U_kn[k,:] *= kJperkcal
        V_kn[k,:] *= nm3perA3
    else:
        print "The file type %s isn't defined!" % (options.filetype)
        sys.exit()

# compute correlation times for the data
# Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
print "Now determining correlation time"
g = numpy.ones(2);
ge = numpy.ones(2);
gv = numpy.ones(2);

if (options.efficiency == None):
    if (type != 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
        for k in range(2):
            ge[k] = timeseries.statisticalInefficiency(U_kn[k,0:N_k[k]],fast=False)
        print "statistical inefficiencies of energy are %.3f and %.3f steps" % (ge[0],ge[1])
    if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
        for k in range(2):
            gv[k] = timeseries.statisticalInefficiency(V_kn[k,0:N_k[k]],fast=False)
        print "statistical inefficiencies of volume are %.3f and %.3f steps" % (gv[0],gv[1])
    for k in range(2):
        g[k] = numpy.max([ge[k],gv[k]])
    print "Using %.3f and %.3f  as the statistical inefficiencies" % (g[0],g[1])
else:
    for k in range(2):
        g[k] = options.efficiency[k]
    print "statistical inefficiencies taken from input options and are %.3f and %.3f steps" % (options.efficiency[0],options.efficiency[1])
if (options.useg == 'subsample'):
    N_k_sampled = numpy.zeros(2)
    tempspace = numpy.zeros(numpy.max(N_k))
    for k in range(2):
        if (type != 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
            indices = timeseries.subsampleCorrelatedData(U_kn[k,0:N_k[k]],g[k])
            tempspace = U_kn[k,indices].copy()
            N_k_sampled[k] = numpy.size(indices) 
            U_kn[k,0:N_k_sampled[k]] = tempspace[0:N_k_sampled[k]]
        if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
            indices = timeseries.subsampleCorrelatedData(V_kn[k,0:N_k[k]],g[k])
            tempspace = V_kn[k,indices].copy()
            N_k_sampled[k] = numpy.size(indices) 
            V_kn[k,0:N_k_sampled[k]] = tempspace[0:N_k_sampled[k]]
        print "data has been subsampled using the statistical inefficiencies"
        g[k] = 1.0
        N_k[k] = N_k_sampled[k]
else:
    print "statistical efficiencies used to scale the statistical uncertained determined from all data"
figname = options.figname
title = options.figname

ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,P_k=P_k,U_kn=U_kn,V_kn=V_kn,title=title,figname=figname,nbins=nbins,
                    reptype='bootstrap',g=g,nboots=nboots,bMaxwell=(type=='kinetic'),bLinearFit=bLinearFit,bNonLinearFit=bNonLinearFit,bMaxLikelihood=bMaxLikelihood,seed=options.seed)
