#!/usr/bin/python

import numpy
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.integrate import quad
from scipy.stats import geom 
import checkdist
import optparse, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--temperature", nargs = 2, dest="T_k", type="float",default=[0.5,1.5],
                  help="low and high temperatures, [default = %default]") 
parser.add_option("-m", "--mu", nargs = 2, dest="mu_k", type="float",default=[0.6,1.4],
                  help="low and high chemical potentials, [default = %default]") 
parser.add_option("-c", "--cuttails", dest="cuttails", type="float",default=0.01,
                  help="fraction of the tails to omit from the analysis to avoid small sample errors in binning [default = %default]")
parser.add_option("-b", "--nboots", dest="nboots", type="int",default=200,
                  help="number of bootstrap samples performed [default = %default]")
parser.add_option("-r", "--nreps", dest="nreps", type="int",default=0,
                  help="number of independent repetitions of the sampling [default = %default]")
parser.add_option("-i", "--nbins", dest="nbins",type = "int", default=30,
                  help="number of bins for bootstrapping [default = %default]") 
parser.add_option("-e", "--energytype", dest="type", default="jointEN",
                  help="the type of energy that is being analyzed [default = %default]")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",default=False,
                  help="more verbosity")
parser.add_option("-N", "--number", nargs = 2, dest="N_k", type="int",default=[50000,50000],
                  help="number of samples from the two states, [default = %default]") 
parser.add_option("-K", "--force", dest="K", type="float",default=1.0,
                  help="spring force constant prefactor [default = %default]")
parser.add_option("-g", "--figureprefix", dest="figname", default='harmonic_mu',
                  help="name prefix for the figure")
parser.add_option("-s", "--seed", dest="seed", type = "int", default=None,
                  help="random seed for generating independet or bootstrap samples")

(options, args) = parser.parse_args()

# For N 1D harmonic oscillator, the potential is given by                                
#   V(x;K) = (K/2) * (x_-x_0)**2                                                                             
# where K denotes the spring constant.                                  
#                                                
# The equilibrium distribution is given analytically by                                      
#   p(x;beta,N) = sqrt[(beta K) / (2 pi)] exp[-beta K (x-x_0)**2 / 2]                                                       
# The dimensionless free energy is therefore                                                
#   q(beta,K) = sqrt(2 pi / beta K)
#   f(beta,K) = - (1/2) * ln[ (2 pi) / (beta K) ]              
#   f_2 - f_1 = -(1/2) (ln [ 1/ beta_2 K_2] - ln [ 1/beta_1 K_1])  
#   f_2 - f_1 = (1/2) ln [beta_2/K_2 / beta_1/K_1]  
#   f_2 - f_1 = (1/2) ln [K_1 T_1/ K_2 T_2]  
#   for T_1 = 0.9, T_2 = 1.1, K_1 = K_2, df = 0.5*ln(1.1/0.9) = 0.10034
#   for T_1 = 0.8, T_2 = 1.2, K_1 = K_2, df = 0.5*ln(1.2/0.8) = 0.2027
#   for T_1 = 0.5, T_2 = 1.5, K_1 = K_2, df = 0.5*ln(1.5/0.5) = 0.5493
#   Now add number of particles.  Use the model of distinguishable noninteracting particles.
#   In this case, Q(beta,N) = Q^N

#   Xi(x,N) = \sum_N Q(beta,N) exp(beta mu N)
#           = \sum_N (2 pi /beta K)^(N/2) exp(beta mu N)
#           = \sum_N sqrt(2 pi /beta K)^N exp(beta mu)^N
#           = \sum_N [sqrt(2 pi /beta K)^N exp(beta mu)]^N
#      for code simplicity, we assume N=1 is the minimum 
#           set x = (sqrt(2 pi / beta K)*exp(beta mu))
#           = ((1-x)^-1)-1 = 1/(1-x) - (1-x)/(1-x) = x/(1-x)
#
#   P(x,N) propto exp(-\beta (E - mu N))
#   P(x,N) exp(-\sum_i \betaK/2 x_i^2 -\beta mu N)/[1-(sqrt(2 pi/ beta K)exp(beta mu))]^(-1)
#   
# If there are two temperatures and the same chemical, we just bin the A_i = E_i - \mu N_i 
# So we want to sample from two temperatures and two mu.

# P_1 / P_2 = exp(B_2 PV_1 - B_2 PV_2 - (B_2 - B_1) E + (B_2 mu_2 - B_2 mu_2) N) 

#   analytical variance = depends on the method. 

title='harmonic oscillators with number'
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
a_k = options.K*numpy.array([1,1])
T_k = numpy.array(options.T_k) #T_k = numpy.array([0.5,1.5])  # temperatures
mu_k = numpy.array(options.mu_k) #T_k = numpy.array([0.5,1.5])  # chemical potentials
beta_k = 1.0/(kB*T_k)  # kB = 1 
N_k = numpy.array(options.N_k) #N_k number of samples

if (options.type == 'helmholtz'):
    analysis_type = 'dbeta-constmu'
    if (mu_k[0] != mu_k[1]):
        print "Chemical potentials are not equal: can't test the ethalpy distribution"
    if (T_k[0] == T_k[1]):
        print "Temperatures are equal: can sometimes result in numerical instability"
elif (options.type == 'number'):
    analysis_type = 'dmu-constB'
    if (T_k[0] != T_k[1]):
        print "Temperatures are not equal: can't test the volume distribution"
    if (mu_k[0] == mu_k[1]):
        print "Chemical potentials are equal: can sometimes result in numerical instability"
elif (options.type == 'jointEN'):
    analysis_type = 'dbeta-dmu'
    if (T_k[0] == T_k[1]):
        print "Temperatures are equal: can sometimes result in numerical instability"
    if (mu_k[0] == mu_k[1]):
        print "Chemical potentials are equal: can sometimes result in numerical instability"
else:
    print "analysis type %s is not defined!" % (options.type)

N_max = numpy.max(N_k)

gtau = 1
genfreq = 1000

K = len(beta_k)
K_k = numpy.array([1,1]) # the spring constants
O_k = numpy.array([0,0]) # the locations
sigma_k = (beta_k * K_k)**(-0.5)
N_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the number of harmonic oscillators in the sample 
U_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the energy of the sample at x_kn[k,n]

#f     = -log [1-(sqrt(2 pi /beta K)exp(-beta mu))]^(-1)
#f     = log [1-(sqrt(2 pi /beta K)exp(-beta mu))]
#df     = log ([1-(sqrt(2 pi /beta[1] K) exp(-beta[1] mu[1]))] / [1-(sqrt(2 pi /beta[0] K)exp(-beta[0] mu[0]))]) 
#           set x = (sqrt(2 pi / beta K)*exp(-beta mu))
#           = ((1-x)^-1)-1 = 1/(1-x) - (1-x)/(1-x) = x/(1-x)
xv = numpy.sqrt(2*numpy.pi/(beta_k*K_k))*numpy.exp(beta_k*mu_k)
f  = -numpy.log(xv/(1-xv))
df = f[1]-f[0]
print "Analytical df = %.8f" % (df)

if (options.nreps > 0):
    print "Now sampling %d sets of data . . . could also take a bit" % (nreps)

reps = []

if (reptype == 'independent'):
    ncount = nreps
elif (reptype == 'bootstrap'):
    ncount = 1

for n in range(ncount):
    if (n%10 == 0 and n>0):
        print "Finished generating %d data sets . . ." % (n)

    # generate independent samples from (x,N).  
    # Pick samples from N, and then from x given N 
    # P(N) propto [Q exp(-\beta mu)]^N \propto [sqrt(2 Pi /beta K)exp(beta mu)]^N  
    # P(X) = a^x 
    # then sample N with the geometric distribution. 

    for k in range(K):

        for n in range(N_k[k]):  # gibbs sampling to get x,N samples:
            
            if (n%10000) == 0:
                print "Set %d: Generated up through sample %d" % (k,n)
            # for the x coordinate

            p = numpy.exp(beta_k[k]*mu_k[k])*numpy.sqrt(2*numpy.pi/beta_k[k]*K_k[k])     
            N_kn[k,n] = scipy.stats.geom.rvs(1-p,size=1)

            # now generate random distances
            x_i = numpy.random.normal(O_k[k], sigma_k[k], N_kn[k,n])

            # compute potential energy of all samples in all potentials 
            U_kn[k,n] = 0.5*(K_k[k] * numpy.sum(x_i**2)) 

    addrep = [U_kn.copy(),N_kn.copy()]    
    reps.append(addrep)

checkdist.ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,mu_k=mu_k,U_kn=U_kn,N_kn=N_kn,kB=1.0,title=title,figname=options.figname,nbins=options.nbins, reptype=reptype, nboots=options.nboots, reps=reps, cuttails=options.cuttails, eunits='kT',seed=options.seed)
