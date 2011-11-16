#!/usr/bin/python

import numpy
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import checkdist
from checkdist import *

import pdb

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
reptype = 'independent'
nreps = 100
nbins = 40
bMaxLikelihood = True
kB = 1.0
K_k = numpy.array([1,1])
T_k = numpy.array([0.5,1.5])
scale = 1.5
T_k = numpy.array([numpy.sqrt(1.0/scale),numpy.sqrt(scale)])
noise = 0.01;  #Random noise

figname='harmonic_' + str(scale) + '_' + str(noise)

beta_k = 1.0/(kB*T_k)  # kB = 1 
O_k = numpy.array([0,0])
N_k = 10000*numpy.array([1,1])
sigma_k = (beta_k * K_k)**(-0.5)
N_max = numpy.max(N_k)

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
        #U_kn[k,:] += numpy.abs(noise*numpy.random.normal(O_k[k],1,[N_k[k]])) # add noise

    addrep = [U_kn.copy()]    
    reps.append(addrep)

ProbabilityAnalysis(N_k,T_k=T_k,U_kn=U_kn,kB=1.0,title=title,figname=figname,nbins=nbins,reptype=reptype,cuttails=0.01, reps=reps,eunits='kT')
