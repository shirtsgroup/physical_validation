#!/usr/bin/python

import numpy
import numpy.random
import matplotlib
import matplotlib.pyplot as plt
import pymbar
import scipy
import scipy.optimize
import scipy.stats
import pdb

#==========================
# HELPER FUNCTIONS
#=========================

def LogLikelihood(x,beta,U_kn,N_k):

    M = numpy.log(N_k[1]/N_k[0])

    D0 = M + beta*x[0] + x[1]*U_kn[0,0:N_k[0]]
    D1 = M + beta*x[0] + x[1]*U_kn[1,1:N_k[0]]
    of = numpy.sum(numpy.log(1 + numpy.exp(D0))) + numpy.sum(numpy.log(1 + numpy.exp(-D1)))
                  
    return of

def MaxLikeParams(beta_k,U_kn,N_k,df=0,analytic_uncertainty=False):

    beta = numpy.average(beta_k)
    trueslope = -(beta_k[1]-beta_k[0])
    optimum = scipy.optimize.fmin(LogLikelihood,[df,trueslope],args=(beta,U_kn,N_k),disp=0);
    optimum[0] *= beta

    if (analytic_uncertainty):
        return optimum, doptimum 
    else:
        return optimum

def printstats(title, kB,T_k,df,ddf,slope,dslope,trueslope):
    print "---------------------------------------"
    print "        %20s        " % (title)
    print "---------------------------------------"
    print ""
    print "---------------------------------------------"
    print "     df = %.5f +/- %.5f " % (df,ddf)
    print "---------------------------------------------"
    print "     Estimated slope       vs.   -(b_2-b_1)"
    print "---------------------------------------------"
    print "%11.6f +/- %11.6f |   %11.6f" % (slope, dslope, trueslope)
    print "---------------------------------------------"

    quant = numpy.abs((slope-trueslope)/dslope)
    print ""
    print "(That's %.2f quantiles from true slope=%5f, FYI.)" % (quant,trueslope)
    if (quant > 3):
        print " (Ouch!)"
    else:
        print ""
    print "---------------------------------------------"
    print " True dT = %7.3f, Eff. dT = %7.3f+/-%.3f" % (T_k[1]-T_k[0],kB*slope*T_k[0]*T_k[1], dslope*kB*T_k[0]*T_k[1])
    print "---------------------------------------------"

def PrintPicture(xaxis,true,y,dy,fit,name,figname,show=False):

    import matplotlib
    import matplotlib.pyplot as plt

    print "Now printing figure"
    plt.clf()
    plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r'$\ln\frac{P_2(E)}{P_1(E)}$')
    plt.errorbar(xaxis,true,fmt='k-',label = r'$-(\beta_2-\beta_1)E$')
    plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to y = aE+b')
    plt.title('E vs. probability ratio for ' + name)
    plt.xlabel('$E$ (kT)')
    plt.ylabel(r'$\ln\frac{P_1(E)}{P_2(E)}$')
    plt.legend(loc='lower right')
    if show:
        plt.show()
    plt.savefig(figname + '.pdf')


def LinFit(beta_k,U_kn,N_k,bins,df=0,analytic_uncertainty=False,bGraph=False,name="",figname='figure.pdf',type='kinetic'):

    K = len(N_k)

    hlist = []
    dhlist = []

    for k in range(0,K):
        hstat = plt.hist(U_kn[k,0:N_k[k]], bins = bins)
        h = (1.0*hstat[0])/N_k[k] 
        hlist.append(h)
        dh = numpy.sqrt(h*(1.0-h)/N_k[k])
        dhlist.append(dh)

    ratio = numpy.log(hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2)

    usedat = numpy.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0/dratio[usedat]

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    
    x = xaxis[usedat]

    Z = numpy.ones([nuse,2],float)
    Z[:,1] = x
    
    # the true line is y = df + -(b_2-b_1) x, where y is ln P_1(E)/P_2(E)

    W = numpy.diag(weights) # weights should be squared?
    M = numpy.dot(W,Z)
    MT = numpy.transpose(M)
    IN = numpy.dot(MT,M)
    WY = numpy.dot(W,y)
    MY = numpy.dot(MT,WY)
    X = numpy.linalg.solve(IN,MY)

    trueslope = -(beta_k[1]-beta_k[0])

    if (bGraph):
        true = df+trueslope*xaxis 
        fit = X[0] + X[1]*xaxis
        print "      X         True     Observed     Error   d(true/obs) sig(true/obs)  Fit   "
        print "---------------------------------------------------------------------------------------"
        for i in range(len(ratio)):
            diff = ratio[i]-true[i]
            sig = numpy.abs(ratio[i]-true[i])/dratio[i]
            print "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f" % (xaxis[i],true[i],ratio[i],dratio[i],diff,sig,fit[i])

        PrintPicture(xaxis,true,ratio,dratio,fit,name,figname)

    if (analytic_uncertainty):
        return X,dX
    else:
        return X

def ProbabilityAnalysis(U_kn,T_k,N_k,kB=0.0083144624,title=None,figname=None,nbins=40,
                        bMaxLikelihood=True,bLinearFit=True,reptype=None,nboots=200,reps=None,cuttails=100):

    K = 2  # should be 2 for now

    beta_k = (1.0/(kB*T_k))

    # determine the bin widths
    maxk = numpy.zeros(K,float)
    mink = numpy.zeros(K,float)

    # cuttails indicates how many we leave out on each tail
    # for now, we choose the range that cuts 100 from the tails of the smallest distribution.

    minN = min(N_k)
    prange = (1.0*cuttails)/minN

    for k in range(K):
        maxk[k] = scipy.stats.scoreatpercentile(U_kn[k,0:N_k[k]],100*(1-prange))
        mink[k] = scipy.stats.scoreatpercentile(U_kn[k,0:N_k[k]],100*(prange))

    binmax = numpy.min(maxk)
    binmin = numpy.max(mink)

    bins = numpy.zeros(nbins+1,float)
    for i in range(nbins+1):
        bins[i] = binmin + (binmax-binmin)*(i/(1.0*nbins))    

    #===================================================================================================
    # Calculate free energies with different methods
    #===================================================================================================    

    trueslope = -(beta_k[1]-beta_k[0])

    print "True slope of log(P_2(E)/P_1(E)) should be %.8f" % trueslope
    print "Now computing -ln(Q_1/Q_2) using BAR"

    w_F = (beta_k[1]-beta_k[0])*U_kn[0,0:N_k[0]];
    w_R = (beta_k[0]-beta_k[1])*U_kn[1,0:N_k[1]];

    (df,ddf) = pymbar.BAR(w_F,w_R,compute_uncertainty=True)

    print "using %.5f for -ln(Q_1/Q_2) computed from BAR" % (df) 
    print "Uncertainty in quantity is %.5f" % (ddf)
    print "Assuming this is negligible compared to sampling error at individual points" 

    print "Now compute the linear fit"

    if (bLinearFit):
        (lindf,linslope) = LinFit(beta_k,U_kn,N_k,bins,df=df,name=title,figname=figname,bGraph=True)

    print "Linear Fit: df = %.8f, slope = %.8f" % (lindf,linslope)

    print "Now compute the maximum likelihood version slope" 

    if (bMaxLikelihood):
        (mldf,mlslope) = MaxLikeParams(beta_k,U_kn,N_k,df=df)

    print "Maximum likelihood: df = %.8f, slope = %.8f" % (mldf,mlslope)

    if (reptype == None):
        return

    if (reptype == 'bootstrap'):
        nreps = nboots
        print "Now bootstrapping (n=%d) for uncertainties . . . could take a bit!" % (nboots)
    elif (reptype == 'independent'):
        print "Now analyzing %d independent samples . . . could also take a bit!" % (nreps)
    else:
        print "Don't understand reptype = %s; quitting" % (reptype)

    linslopes = numpy.zeros([nreps],float)
    lindfs = numpy.zeros([nreps],float)
    mlslopes = numpy.zeros([nreps],float)
    mldfs = numpy.zeros([nreps],float)
    if (reptype == 'bootstrap'):    
        Ur_kn = numpy.zeros([K,numpy.max(N_k)],float)

    for n in range(nreps):
        if (n%10 == 0):
            print "Finished %d samples . . ." % (n)
        
        if (reptype == 'bootstrap'):    
            for k in range(K):
                rindex = numpy.random.randint(0,high=N_k[k],size=N_k[k]);  # bootstrap it 
                Ur_kn[k,0:N_k[k]] = U_kn[k,rindex]

        if (reptype == 'independent'):
            Ur_kn = reps[n] 

        if (bLinearFit):    
            (lindfs[n],linslopes[n]) = LinFit(beta_k,Ur_kn,N_k,bins) 

        if (bMaxLikelihood):
            (mldfs[n],mlslopes[n]) = MaxLikeParams(beta_k,Ur_kn,N_k,df=df)
            
    if (bLinearFit):
        df = numpy.average(lindfs)
        ddf  = numpy.std(lindfs)
        slope = numpy.average(linslopes)
        dslope = numpy.std(linslopes)

    printstats('Linear Fit Analysis',kB,T_k,df,ddf,slope,dslope,trueslope)

    if (bMaxLikelihood):
        df = numpy.average(mldfs)
        ddf  = numpy.std(mldfs)
        slope = numpy.average(mlslopes)
        dslope = numpy.std(mlslopes)

    printstats('Maximum Likelihood Analysis',kB,T_k,df,ddf,slope,dslope,trueslope)        

    return
    
