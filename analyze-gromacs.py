
#!/usr/bin/python

# Example illustrating the use of MBAR for computing the hydration free energy of OPLS 3-methylindole
# in TIP3P water through alchemical free energy simulations.

#===================================================================================================
# IMPORTS
#===================================================================================================

import numpy
import timeseries
import pdb
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
parser.add_option("-u", "--useefficiency", dest="useg", type = "string", default='scale',
                  help= "calculate the efficiency by scaling the uncertainty, or by subsampling the input data")

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
pdb.set_trace()
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
kB = 1.3806488*6.0221415/1000.0  # Boltzmann's constant (kJ/mol/K)  
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
# allocate space
U_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n
V_kn = numpy.zeros([K,N_max], dtype=numpy.float64) # U_kn[k,n] is the energy of the sample k,n

for k,T in enumerate(T_k):

    # Read contents of file into memory.
    print "Reading %s..." % filenames[k]
    infile = open(filenames[k], 'r')
    lines = infile.readlines()
    infile.close()

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
           time = float(elements[0])
           if (type != 'volume'):
               energy = float(elements[ecol])
               U_kn[k,N_k[k]] = energy   
           if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
               volume = float(elements[vcol])
               V_kn[k,N_k[k]] = volume   
           N_k[k] += 1 

# compute correlation times for the data
# Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
print "Now determining correlation time"
g = numpy.ones(2);
ge = numpy.ones(2);
gv = numpy.ones(2);

if (options.efficiency == None):
    for k in range(2):
        if (type != 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
            ge[k] = timeseries.statisticalInefficiency(U_kn[k,0:N_k[k]],fast=True)
        if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
            ge[k] = timeseries.statisticalInefficiency(V_kn[k,0:N_k[k]])
        g[k] = numpy.max(ge[k],gv[k])
    print "statistical inefficiencies are %.3f and %.3f steps" % (g[0],g[1])
else:
    for k in range(2):
        g[k] = options.efficiency[k]
    print "statistical inefficiencies taken from input options and are %.3f and %.3f steps" % (options.efficiency[0],options.efficiency[1])
if (options.useg == 'subsample'):
    tempspace = numpy.zeros(numpy.max(N_k))
    if (type != 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
        indices = timeseries.subsampleCorrelatedData(U_kn[k,0:N_k[k]],g[k])
        tempspace = U_kn[k,indices].copy()
        N_k[k] = numpy.size(indices) 
        U_kn[k,0:N_k[k]] = tempspace[0:N_k[k]]
    if (type == 'volume') or (type == 'enthalpy') or (type == 'jointEV'):
        indices = timeseries.subsampleCorrelatedData(V_kn[k,0:N_k[k]],g[k])
        tempspace = V_kn[k,indices].copy()
        N_k[k] = numpy.size(indices) 
        V_kn[k,0:N_k[k]] = tempspace[0:N_k[k]]
    print "data has been subsampled using these statistical inefficiencies"
    g[0] = g[1] = 1.0

figname = options.figname
title = options.figname

ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,P_k=P_k,U_kn=U_kn,V_kn=V_kn,title=title,figname=figname,nbins=nbins,
                    reptype='bootstrap',g=g,nboots=nboots,bMaxwell=(type=='kinetic'),seed=options.seed)
