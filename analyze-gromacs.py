#!/usr/bin/python

# Example illustrating the use of MBAR for computing the hydration free energy of OPLS 3-methylindole
# in TIP3P water through alchemical free energy simulations.

#===================================================================================================
# IMPORTS
#===================================================================================================

import numpy
import pdb
import checkdist
from checkdist import *
import optparse, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--files", dest="datafiles",nargs = 2,
                  help="the two files of different temperature for analysis")
parser.add_option("-d", "--directory", dest="datafile_directory",default ='./',
                  help="the directory the data files are is in")
parser.add_option("-k", "--nolikelihood", dest="bMaxLikelihood", action="store_false",default=True,
                  help="Don't run maximum likelihood analysis [default = run this analysis]") 
parser.add_option("-n", "--nolinearfit", dest="bLinearFit", action="store_false",default=True,
                  help="Don't run linear fit analysis [default = run this analysis]") 
parser.add_option("-t", "--temperature", nargs = 2, dest="T_k", type="float",default=[295,305],
                  help="low and high temperatures, [default = %default]") 
parser.add_option("-e", "--energytype", dest="type", default="total",
                  help="the type of energy that is being analyzed [default = %default]")
parser.add_option("-b", "--nboot", dest="nboots", type="int",default=200,
                  help="number of bootstrap samples performed [default = %default]")
parser.add_option("-i", "--nbins", dest="nbins",type = "int", default=30,
                  help="number of bins for bootstrapping [default = %default]") 
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",default=False,
                  help="more verbosity")
parser.add_option("-p", "--picturename", dest="figname", default='figure.pdf',
                  help="name for the figure")

(options, args) = parser.parse_args()

type = options.type
if (type != 'kinetic') and (type != 'potential') and (type != 'total'):
    print "type of energy %s isn't defined!" % (type)
    sys.exit()

#===================================================================================================
# CONSTANTS
#===================================================================================================

verbose = options.verbose
T_k = numpy.array(options.T_k) #T_k = numpy.array([132.915071475571,137.138128524429])  # temperatures
names = ['down','up']
type = options.type # 'potential', 'kinetic' or 'total'
nboots = options.nboots
nbins = options.nbins
bMaxLikelihood = options.bMaxLikelihood
bLinearFit = options.bLinearFit
figname = options.figname

if (verbose):
    print "verbosity is %s" % (str(verbose))
    print "Energy type is %s" % (type)
    print "\'%s\' temperature is %f" % (names[0],T_k[0])
    print "\'%s\' temperature is %f" % (names[1],T_k[1])
    print "Number of bootstraps is %d" % (nboots)
    print "Number of bins (only used for linear fit) is %d" % (nbins)
    if (bMaxLikelihood):
        print "Generating maximum likelihood statistics"
    else:
        print "Not generating maximum likelihood statistics"
    if (bLinearFit):
        print "Generating linear fit statistics"
    else:
        print "Not generating linear fit statistics"
    print "Figures will be names %s" % (figname)    

# Shouldn't need to modify below this for standard usage 
# ------------------------
K = 2;
kB = 1.3806488*6.0221415/1000.0  # Boltzmann's constant (kJ/mol/K)  
N_k = numpy.zeros([K],int) # number of samples at each state

# check just size of all files
N_size = numpy.zeros(K,int) 
filenames = []
for k,T in enumerate(T_k):
    filename = options.datafile_directory + options.datafiles[k]
    filenames.append(filename)
    print "checking size of \'%s\' temperature file %s..." % (names[k],filenames[k])    
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    N_size[k] = len(lines)

N_max = numpy.max(N_size)
# allocate space
U_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the energy of the sample at x_kn[k,n]

for k,T in enumerate(T_k):

    # Read contents of file into memory.
    print "Reading %s..." % filenames[k]
    infile = open(filenames[k], 'r')
    lines = infile.readlines()
    infile.close()

    match = False
    for line in lines:
        # Split line into elements.
        if (line[0:3] == '@ s'):
            elements = line.split()
            whichcol = int((elements[1])[1:])+1   # figure out which column it is
            if (type == 'potential'):
                if (elements[3] == "\"Potential\""):
                    col = whichcol
                    match = True
            if (type == 'total'):
                if (elements[3] == "\"Total"):
                    comp = elements[3] + ' ' + elements[4]
                    if (comp == "\"Total Energy\""):
                        col = whichcol
                        match = True
            if (type == 'kinetic'):
                if (elements[3] == "\"Kinetic"):
                    comp = elements[3] + ' ' + elements[4]
                    if (comp == "\"Kinetic En.\""):
                        col = whichcol
                        match = True

        if ((line[0] != '#') and (line[0] != '@')):
                
           elements = line.split()
           # what is the time of the sample
           time = float(elements[0])
           energy = float(elements[col])
           U_kn[k,N_k[k]] = energy   
           N_k[k] += 1 

figname = options.figname
title = options.figname
ProbabilityAnalysis(U_kn,T_k,N_k,title=title,figname=figname,nbins=nbins,
                    bMaxLikelihood=True,bLinearFit=True,reptype='bootstrap',nboots=nboots)
