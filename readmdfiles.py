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

onlyE = ['potential', 'kinetic', 'total']
requireV = ['enthalpy', 'volume', 'jointEV'] 
requireN = ['helmholtz', 'number', 'jointEN']
alltypes = onlyE + requireV + requireN

def read_flatfile(lines,type,N_max):

    # assumes kJ/mol energies, nm^3 volumes 
    # allocate space

    U_n = numpy.zeros([N_max], dtype=numpy.float64) # U_n[k,n] is the energy of the sample n
    V_n = numpy.zeros([N_max], dtype=numpy.float64) # V_n[k,n] is the volume of the sample n
    N_n = numpy.zeros([N_max], dtype=numpy.float64) # N_n[k,n] is the number of particles of the sample n
    N = 0

    # we assume energy is first, then volume, then n
    for line in lines:
        if (line[0] != '#'):   # in flat file format, anything that starts with a hash is an ignored comment
        
            elements = line.split()
            numcol = len(elements)
            if (numcol == 0):
                print "Error: No data for data point %d" % (N) 
                sys.exit()
            elif (numcol == 1):     
                if (type in onlyE):
                    U_n[N] = float(elements[0])                    
                elif (type == 'volume'):
                    V_n[N] = float(elements[0])                    
                elif (type == 'number'):
                    N_n[N] = float(elements[0])                    
                else:
                    print "Error: asking for test requiring multiple variables (%s) but only provided one column of data" % (type)
                    sys.exit()
            elif (numcol == 2):
                if type in requireV and type != 'volume':
                    U_n[N] = float(elements[0])                    
                    V_n[N] = float(elements[1])                    
                elif type in requireN and type != 'number':
                    U_n[N] = float(elements[0])                    
                    N_n[N] = float(elements[1])                    
                else:
                    print "Error: asking for test (%s) incompatible with two columns of data" % (type)
            elif (numcol > 2):
                print "Error: there is no test that required the provided %d columns of data" % (numcol)
                sys.exit()
            N += 1 

    return U_n,V_n,N_n,N        

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

def getefficiency(N_k,U_kn,V_kn,N_kn,type):

    K = len(N_k)
    g = numpy.ones(K)
    ge = numpy.ones(K);
    gv = numpy.ones(K);
    gn = numpy.ones(K);

    if (type != 'volume') and (type != 'number'):
        for k in range(K):
            ge[k] = timeseries.statisticalInefficiency(U_kn[k,0:N_k[k]],fast=False)
        print "Calculating ["
        for k in range(K):
            print " %.3f " % (ge[k])
        print "] as the statistical inefficiencies of the energy"
    if type in requireV:
        for k in range(K):
            gv[k] = timeseries.statisticalInefficiency(V_kn[k,0:N_k[k]],fast=False)
        print "Calculating ["
        for k in range(K):
            print " %.3f " % (gv[k])
        print "] as the statistical inefficiencies of the volume"
    if type in requireN:
        for k in range(K):
            gn[k] = timeseries.statisticalInefficiency(N_kn[k,0:N_k[k]],fast=False)
        print "Calculating ["
        for k in range(K):
            print " %.3f " % (gn[k])
        print "] as the statistical inefficiencies of the particle number"

    for k in range(K):
        g[k] = numpy.max([ge[k],gv[k],gn[k]])
    print "Using ["
    for k in range(K):
        print " %.3f " % (g[k])
    print "] as the statistical inefficiencies"
    return g

def subsample(N_k,U_kn,V_kn,N_kn,g,type):

    K = len(N_k)
    N_k_sampled = numpy.zeros(K)
    tempspace = numpy.zeros(numpy.max(N_k))
    for k in range(K):
        if (type != 'volume') and (type != 'number'): 
            indices = timeseries.subsampleCorrelatedData(U_kn[k,0:N_k[k]],g[k])
            tempspace = U_kn[k,indices].copy()
            N_k_sampled[k] = numpy.size(indices) 
            U_kn[k,0:N_k_sampled[k]] = tempspace[0:N_k_sampled[k]]
        if (type in requireV):
            indices = timeseries.subsampleCorrelatedData(V_kn[k,0:N_k[k]],g[k])
            tempspace = V_kn[k,indices].copy()
            N_k_sampled[k] = numpy.size(indices) 
            V_kn[k,0:N_k_sampled[k]] = tempspace[0:N_k_sampled[k]]
        if (type in requireN):
            indices = timeseries.subsampleCorrelatedData(N_kn[k,0:N_k[k]],g[k])
            tempspace = N_kn[k,indices].copy()
            N_k_sampled[k] = numpy.size(indices) 
            N_kn[k,0:N_k_sampled[k]] = tempspace[0:N_k_sampled[k]]
        print "data has been subsampled using the statistical inefficiencies"
        g[k] = 1.0
        N_k[k] = N_k_sampled[k]
