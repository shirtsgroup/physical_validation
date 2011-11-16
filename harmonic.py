#!/usr/bin/python

import numpy
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import checkdist
from checkdist import *

import optparse, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--temperature", nargs = 2, dest="T_k", type="float",default=[0.5,1.5],
                  help="low and high temperatures, [default = %default]") 
parser.add_option("-c", "--cuttails", dest="cuttails", type="float",default=0.001,
                  help="fraction of the tails to omit from the analysis to avoid small sample errors in binning [default = %default]")
parser.add_option("-K", "--force", dest="K", type="float",default=1.0,
                  help="spring force constant[default = %default]")
parser.add_option("-b", "--nboot", dest="nboots", type="int",default=0,
                  help="number of bootstrap samples performed [default = %default]")
parser.add_option("-r", "--nreps", dest="nreps", type="int",default=200,
                  help="number of independent repetitions of the sampling [default = %default]")
parser.add_option("-i", "--nbins", dest="nbins",type = "int", default=30,
                  help="number of bins for bootstrapping [default = %default]") 
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",default=False,
                  help="more verbosity")
parser.add_option("-N", "--number", nargs = 2, dest="N_k", type="int",default=[100000,100000],
                  help="number of samples from the two states, [default = %default]") 
parser.add_option("-g", "--figureprefix", dest="figname", default='harmonic',
                  help="name prefix for the figure")
parser.add_option("-o", "--noise", dest="noise", type = 'float', default=0.0,
                  help="random noise perturbing the potential")
parser.add_option("-s", "--seed", dest="seed", type = 'int', default=None,
                  help="random seed for generating independent or bootstrap samples")


(options, args) = parser.parse_args()

# For a 1D harmonic oscillator, the potential is given by                                
#   V(x;K) = (K/2) * (x-x_0)**2                                                                             
# where K denotes the spring constant.                                  
#                                                                                                   
# The equilibrium distribution is given analytically by                                      
#   p(x;beta,K) = sqrt[(beta K) / (2 pi)] exp[-beta K (x-x_0)**2 / 2]                                                       
# The dimensionless free energy is therefore                                                
#   f(beta,K) = - (1/2) * ln[ (2 pi) / (beta K) ]              
#   f_2 - f_1 = -(1/2) (ln [ 1/ beta_2 K_2] - ln [ 1/beta_1 K_1])  
#   f_2 - f_1 = (1/2) ln [beta_2/K_2 / beta_1/K_1]  
#   f_2 - f_1 = (1/2) ln [K_1 T_1/ K_2 T_2]  
#   for T_1 = 0.9, T_2 = 1.1, K_1 = K_2, df = 0.5*ln(1.1/0.9) = 0.10034
#   for T_1 = 0.8, T_2 = 1.2, K_1 = K_2, df = 0.5*ln(1.2/0.8) = 0.2027
#   for T_1 = 0.5, T_2 = 1.5, K_1 = K_2, df = 0.5*ln(1.5/0.5) = 0.5493
#   analytical variance = depends on the method.  So won't try 

title='harmonic oscillators'
if (options.nreps > 0):
    reptype = 'independent'
    nreps = options.nreps
if (options.nboots > 0):
    reptype = 'bootstrap'
    nreps = options.nboots
if (options.nboots > 0 and options.nreps > 0):
    print "Can't do both bootstrap sampling and independence sampling: defaulting to bootstrap sampling"

if (options.seed):
    numpy.random.seed(options.seed) # setting the seed for independent sampling 
    print "setting random number seed for generating samples as %d" % (options.seed)

kB = 1.0
K_k = options.K*numpy.array([1,1])
T_k = numpy.array(options.T_k) #T_k = numpy.array([0.5,1.5])  # temperatures
noise = options.noise  #Random noise
beta_k = 1.0/(kB*T_k)  # kB = 1 
O_k = numpy.array([0,0])
N_k = numpy.array(options.N_k) #N_k number of samples
sigma_k = (beta_k * K_k)**(-0.5)
N_max = numpy.max(N_k)

if (T_k[0] == T_k[1]):
    print "Temperatures are equal: can sometimes result in numerical instability"

K = len(beta_k)
x_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the coordinate x of independent snapshot n of simulation k   
U_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the energy of the sample at x_kn[k,n]

df = 0.5*numpy.log(beta_k[1]/beta_k[0]) # analytical result
print "Analytical df = %.8f" % (df)

print "Now sampling %d sets of data . . . could also take a bit" % (nreps)

reps = []

for n in range(nreps):
    if (n%10 == 0):
        print "Finished generating %d sets . . ." % (n)

    for k in range(0,K):
        # generate N_k[k] independent samples with spring constant K_k[k]                   
        x_kn[k,0:N_k[k]] = numpy.random.normal(O_k[k], sigma_k[k], N_k[k])
        # compute potential energy of all samples in all potentials 
        U_kn[k,:] = (K_k[k]*(1+noise)/2.0) * (x_kn[k,:]-O_k[k])**2 
        if (noise > 0):
            U_kn[k,:] += numpy.abs(noise*numpy.random.normal(O_k[k],1,[N_k[k]])) # add noise to the potential

    addrep = [U_kn.copy()]    
    reps.append(addrep)

ProbabilityAnalysis(N_k,T_k=T_k,U_kn=U_kn,kB=1.0,title=title,figname=options.figname,nbins=options.nbins,reptype=reptype,cuttails=options.cuttails, reps=reps,eunits='kT',seed=options.seed)
# OK to pass the same seed, because it will be used for completely different things 
