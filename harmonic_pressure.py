#!/usr/bin/python

import numpy
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.integrate import quad
import checkdist
from checkdist import *

# For a 1D harmonic oscillator, the potential is given by                                
#   V(x;K) = (K/2) * (x-x_0)**2                                                                             
# where K denotes the spring constant.                                  
#                                                                                                   
# The equilibrium distribution is given analytically by                                      
#   p(x;beta,K) = sqrt[(beta K) / (2 pi)] exp[-beta K (x-x_0)**2 / 2]                                                       
# The dimensionless free energy is therefore                                                
#   q(beta,K) = sqrt(2 pi / beta K)
#   f(beta,K) = - (1/2) * ln[ (2 pi) / (beta K) ]              
#   f_2 - f_1 = -(1/2) (ln [ 1/ beta_2 K_2] - ln [ 1/beta_1 K_1])  
#   f_2 - f_1 = (1/2) ln [beta_2/K_2 / beta_1/K_1]  
#   f_2 - f_1 = (1/2) ln [K_1 T_1/ K_2 T_2]  
#   for T_1 = 0.9, T_2 = 1.1, K_1 = K_2, df = 0.5*ln(1.1/0.9) = 0.10034
#   for T_1 = 0.8, T_2 = 1.2, K_1 = K_2, df = 0.5*ln(1.2/0.8) = 0.2027
#   for T_1 = 0.5, T_2 = 1.5, K_1 = K_2, df = 0.5*ln(1.5/0.5) = 0.5493
#   Now add pressure.  Use the model K = (a/V)^2, Then:
#   Xi = \int dV Q(beta,V) exp(-beta PV)
#      = \int dV sqrt(2 pi /beta K) exp( - beta PV)
#      = \int dV sqrt(2 pi V^2 /beta a^2) exp( - beta PV)
#      = sqrt(2 pi/beta a^2) \int dV V exp( - beta PV)
#      = (2 pi/beta)^{1/2} a^{-1} (beta P)^{-2} 
#      = (2 pi)^{1/2} a^{-1} beta^{-5/2} P^{-2} 
#
#   P(x,V) propto exp(-\beta E + -\beta PV)
#          propto exp(-\frac{\beta a^2}{2V^2 + -\beta PV)
#          propto exp(-\beta P[ \frac{a^2}{2PV^2} + V])

# If there are two temperatures and the same pressure, we just bin the enthalpy
# So we want to sample from two temperatures and two pressures.

# P_1 / P_2 = exp(B_2 G_1 - B_2 G_2 - (B_2 - B_1) E - (B_2 P_2 - B_2 P_2) V) 


#   analytical variance = depends on the method. 

title='harmonic oscillators'
figname='harmonic'
#reptype = 'independent'
reptype = 'bootstrap'
nreps = 1
nbins = 30
bMaxLikelihood = True
kB = 1.0
a_k = numpy.array([1,1])

analysis_type = 'dpressure-constB'
#analysis_type = 'dbeta-dpressure'
#analysis_type = 'dbeta-constP'
T_k = numpy.array([1.0,1.0])
#T_k = numpy.array([0.8,1.2])
#T_k = numpy.array([0.5,1.5])
beta_k = 1.0/(kB*T_k)  # kB = 1 
P_k = 1.0*numpy.array([0.7,1.3]) 
#P_k = 1.0*numpy.array([1.0,1.0]) 
N_k = 100000*numpy.array([1,1])
N_max = numpy.max(N_k)

gtau = 1
genfreq = 1000

K = len(beta_k)
x_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the coordinate x of independent snapshot n of simulation k   
V_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the volume of the sample 
U_kn = numpy.zeros([K,N_max], float) # x_kn[k,n] is the energy of the sample at x_kn[k,n]

noise = 0.0;  #Random noise

#df     = -log [(2 pi/beta)^{1/2} a^{-1} beta^{-5/2} P^{-2} / (2 pi/beta)^{1/2} a^{-1} beta^{-5/2} P^{-2}]
#df     = -log [(beta)^{-5/2} (P)^{-2} / (beta)^{-5/2} (P)^{-2}]
df      = 2.5*numpy.log(beta_k[1]/beta_k[0]) + 2.0*numpy.log(P_k[1]/P_k[0]) # analytical result
print "Analytical df = %.8f" % (df)

print "Now sampling %d sets of data . . . could also take a bit" % (nreps)

reps = []

for n in range(nreps):
    if (n%10 == 0):
        print "Finished generating %d sets . . ." % (n)

    # generate N_k[k] independent samples from (x,V)                   
    V = 1; # start with arbitrary volume     

    for k in range(K):

        pfac = (beta_k[k]/2)*(a_k[k])**2
        scale = 1/(beta_k[k]*P_k[k])

        #compute the  envelope constant for this choice of pfac = 
        # k2/k1 = \int (-bV) / int (-aV^(-2) -bV).  The first integral is 1/b, the second
        # must be computed numerically.
        #k2 = scale;
        #k1,dk1 = quad(lambda V: numpy.exp(-pfac*V**(-2)-V/scale),0,numpy.Inf) 
        # this is the scaling factor such that the new distribution is always less than the 
        # one we know exactly
        #M = k2/k1
                  
        for n in range(N_k[k]):  # gibbs sampling to get x,V samples:
            
            if (n%1000) == 0:
                print "Generated up through sample %d" % n
            # for the x coordinate
            gx = numpy.random.normal(0, 1, gtau)

            # for the V coordinate
            ex = numpy.random.exponential((beta_k[k]*P_k[k])**(-1),genfreq)
            ux = numpy.random.random_sample(genfreq)
            j = 0
            for i in range(gtau):  # have to figure this out by experiment . . . 
                
                Kf = (a_k[k]/V)**2
                sigma = (beta_k[k] * Kf)**(-0.5)
                x = sigma*gx[i] 
                while(1):
                    if (j%genfreq==0):
                        # get some more random variates for the V coordinate
                        ex = numpy.random.exponential(scale,genfreq)
                        ux = numpy.random.random_sample(genfreq)
                    j += 1
                    jm = (j-1)%genfreq
                    V = ex[jm]
                    rejection_factor = numpy.exp(-pfac*(x/V)**2)
                    if (ux[jm] < rejection_factor):
                        break;

            x_kn[k,n] = x
            V_kn[k,n] = V
        # generate noise in the energy
        # compute potential energy of all samples in all potentials 
        U_kn[k,:] = 0.5*(a_k[k]/V_kn[k,:])**2 * (x_kn[k,:])**2 

    addrep = [U_kn.copy(),V_kn.copy()]    
    reps.append(addrep)

# compute correlation times for the data
# Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
#import timeseries

# Commenting out, since they seem to be uncorrelated when examined individually.  Still need to look 
# at cross correlation . . . . 
#print "Now determining correlation time"
#g = numpy.ones(2);
#ge = numpy.ones(2);
#gv = numpy.ones(2);
#for k in range(2):
#    ge[k] = timeseries.statisticalInefficiency(V_kn[k,0:N_k[k]])
#    ge[k] = timeseries.statisticalInefficiency(U_kn[k,0:N_k[k]])
#    g[k] = numpy.max(ge[k],gv[k])
#
#print "correlation times are %.3f and %.3f steps" % (g[0],g[1])

#ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,P_k=P_k,U_kn=U_kn,V_kn=V_kn,kB=1.0,title=title,figname=figname,nbins=nbins,reptype=reptype,cuttails=0.02, reps=reps,eunits='kT', vunits="kT", punits="kT")

ProbabilityAnalysis(N_k,type=analysis_type,T_k=T_k,P_k=P_k,U_kn=U_kn,V_kn=V_kn,kB=1.0,title=title,figname=figname,nbins=nbins,reptype='bootstrap',nboots=100,cuttails=0.02,eunits='kT', vunits="kT", punits="kT")
