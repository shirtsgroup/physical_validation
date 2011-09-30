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

def LogLikelihood(x,N_k,beta,U_kn):

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log(N1/N0)

    U0 = U_kn[0,0:N0]
    U1 = U_kn[1,0:N1]

    D0 = M + beta*x[0] + x[1]*U0
    D1 = M + beta*x[0] + x[1]*U1
    E0 = 1 + numpy.exp(D0)
    E1 = 1 + numpy.exp(-D1)

    # this is the negative of the log likelihood, since we want to maximize it using fmin

    of = ((numpy.sum(numpy.log(E0)) + numpy.sum(numpy.log(E1))))/N
    #print of

    return of

def dLogLikelihood(x,N_k,beta,U_kn):

    """
    Derivative with respect to the parameters, to aid the minimization.
    
    """

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log(N1/N0)

    U0 = U_kn[0,0:N0]
    U1 = U_kn[1,0:N1]

    D0 = M + beta*x[0] + U0*x[1]
    D1 = M + beta*x[0] + U1*x[1]

    g = numpy.zeros(2,dtype=numpy.float64)
    E0 = 1/(1 + numpy.exp(-D0))
    E1 = 1/(1 + numpy.exp(D1))

    #this is the gradient of -log likelihood
    g[0] = (beta/N)*(numpy.sum(E0) - numpy.sum(E1))
    g[1] = (1.0/N)*(numpy.sum(U0*E0) - numpy.sum(U1*E1))
                  
    return g

def d2LogLikelihood(x,N_k,beta,U_kn):

    """

    should be I^-1 - (n_f/N + n_r/N), where I is the Fischer information, which is the second
    derivative matrix, evaluated at the maximum likelihood.

    if D = M + beta*x[0] + x[1]*U
    I = \sum_{i=1}^N [[-beta^2/S,-beta*U/S],[-beta*U/S,-U^2/S]] where S = [(1+exp(-D))*(1+exp(D))]

    """

    M = numpy.log(N_k[1]/N_k[0])

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    U = numpy.zeros(N,dtype=numpy.float64)
    U[0:N0] = U_kn[0,0:N0]
    U[N0:N] = U_kn[1,0:N1]

    D = M + beta*x[0] + U*x[1]
    
    E = (1 + numpy.exp(-D)) * (1 + numpy.exp(D))

    hf = numpy.zeros([2,2,N],dtype=numpy.float64)
    
    hf[0,0,:] = beta**2 * numpy.ones(N,dtype=numpy.float64)
    hf[1,0,:] = beta*U
    hf[0,1,:] = beta*U
    hf[1,1,:] = U*U

    # this is the hessian of the minimum function (not the max)
    h = -numpy.sum(hf/E,axis=2)/N
                            
    return h

def solveminlike(x, N_k, beta, U_kn, tol = 1e-10, maxiter=20):
    
    converge = False
    itol = 1e-2
    rtol = 1e-2
    lasttol = 100000;

    for i in range(maxiter):
        lastx = x
        gx = numpy.transpose(dLogLikelihood(x,N_k,beta,U_kn))
        nh = d2LogLikelihood(x,N_k,beta,U_kn)
        dx = numpy.linalg.solve(nh,gx)
        x += dx  # adding, not subtracting because of the handling of negatives 
        rx = dx/x
        checktol = numpy.sqrt(numpy.dot(dx,dx))
        checkrtol = numpy.sqrt(numpy.dot(rx,rx))    
        if (checkrtol < tol):
            break
            converge = True
        if (checkrtol > 1.0) and (checktol > lasttol):  # we are possibly diverging. Switch to cg for a bit.
            x = scipy.optimize.fmin_cg(LogLikelihood,lastx,fprime=dLogLikelihood,gtol=itol,args=[N_k,beta,U_kn],disp=1)
            itol *= rtol
        lasttol = checktol

    if (i == maxiter) and (converge == False):
        print "Too many iterations, convergence failing"

    return x    

def MaxLikeUncertain(x,N_k,beta,U_kn,Uave):

    d = numpy.zeros(2,float)

    # multiply back by N, since we were dealing with smaller numbers for numerical robustness.
    fi = -(N_k[0] + N_k[1])*d2LogLikelihood(x,N_k,beta,U_kn-Uave)

    d2 = numpy.linalg.inv(fi)

    # We have a fit to the line y = m(x-Uave) + b, so to get the uncertainty in the free energy back, we need 
    # to add M*Uave back to f.  The uncertainty in cov(b + m*Uave,b+m*Uave) = var(b) + Uave**2*var(m) + Uave*cov(v,m)

    d[0] = numpy.sqrt(beta*2*d2[0,0] + Uave**2*d2[1,1] - 2*Uave*d2[0,1])
    d[1] = numpy.sqrt(d2[1,1])

    return d

def MaxLikeParams(N_k,beta_k,U_kn,df=0,analytic_uncertainty=False,g=1):

    optimum = numpy.zeros(2,float)
    beta = numpy.average(beta_k)
    trueslope = -(beta_k[1]-beta_k[0])

    # for numerical stability, we need to translate the curve
    Uave = (numpy.sum(U_kn[0,0:N_k[0]]) + numpy.sum(U_kn[1,0:N_k[1]]))/numpy.sum(N_k)

    Umod = U_kn - Uave

    ofit = solveminlike([(df+Uave*trueslope)/beta,trueslope],N_k,beta,Umod,tol=1e-10)
    optimum[0] = beta*(ofit[0])-ofit[1]*Uave # go from free energies to the fit constant
    optimum[1] = ofit[1]

    if (analytic_uncertainty):
        # incorporate g in an average way.
        doptimum = MaxLikeUncertain(ofit,N_k,beta,U_kn,Uave)*numpy.sqrt(1.0/numpy.average(g))
        return optimum[0], doptimum[0], optimum[1], doptimum[1]
    else:
        return optimum[0], optimum[1]

def PrintStats(title,dfs,slopes,kB,T_k,trueslope,ddf='N/A',dslope='N/A'):

    df = numpy.average(dfs) # true even if there is only one
    if (numpy.size(dfs) > 1):
        ddf  = numpy.std(dfs)
        
    slope = numpy.average(slopes) # true even if there is only one
    if (numpy.size(slopes) > 1):    
        dslope = numpy.std(slopes)

    print ""
    print "---------------------------------------------"
    print "     %20s        " % (title)
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

def PrintPicture(xaxis,true,y,dy,fit,name,figname,fittype,show=False):

    import matplotlib
    import matplotlib.pyplot as plt

    print "Now printing figure %s" % (figname)
    plt.clf()
    plt.xlabel('$E$ (kT)')
    if (fittype == 'linear'):
        plt.title('E vs. log probability ratio for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r'$\ln\frac{P_2(E)}{P_1(E)}$')
        plt.errorbar(xaxis,true,fmt='k-',label = r'$-(\beta_2-\beta_1)E$')
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to $y = b+aB$')
        plt.ylabel(r'$\ln\frac{P_1(E)}{P_2(E)}$')
    elif (fittype == 'nonlinear'):
        plt.title('E vs. probability ratio for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r'$\frac{P_2(E)}{P_1(E)}$')
        plt.errorbar(xaxis,true,fmt='k-',label = r'$\exp([\beta_2 F_2- \beta_1 F_1] -[\beta_2-\beta_1]E)$')
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to $y = \exp(b+aE)$')
        plt.ylabel(r'$\frac{P_1(E)}{P_2(E)}$')
    elif (fittype == 'maxwell'):
        plt.title('E vs. probability for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r'$P(E_{\mathrm{kin}})$')
        if (true != None):  # sometimes, true will be none.
            plt.errorbar(xaxis,true,fmt='k-',label = 'Fit to Analytical')
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to Normal')
        plt.ylabel(r'$P(E_{\mathrm{kin}})$')
    else:
        print "I'm crying foul!  Illegal chart type!"

    plt.legend(loc='upper left')
    if show:
        plt.show()
    plt.savefig(figname + '.pdf')

def LinFit(bins,N_k,beta_k,U_kn,df=0,analytic_uncertainty=False,bGraph=False,name="",figname='lin_figure',type='kinetic',g=[1,1]):

    hlist = []
    dhlist = []

    for k in range(0,2):
        hstat = plt.hist(U_kn[k,0:N_k[k]], bins = bins)
        h = (1.0*hstat[0])/N_k[k] 
        hlist.append(h)
        dh = numpy.sqrt(g[k]*h*(1.0-h)/N_k[k])
        dhlist.append(dh)

    ratio = numpy.log(hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2)

    usedat = numpy.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0/dratio[usedat]

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    
    x = xaxis[usedat]

    
    X = numpy.ones([nuse,2],float)
    X[:,1] = x
    
    # the true line is y = df + -(b_2-b_1) x, where y is ln P_1(E)/P_2(E)

    w = numpy.diag(weights) 
    WX = numpy.dot(w,X)
    WY = numpy.dot(w,y)
    WXT = numpy.transpose(WX)
    Z = numpy.dot(WXT,WX)
    WXY = numpy.dot(WXT,WY)

    a = numpy.linalg.solve(Z,WXY)
    da_matrix = numpy.transpose(numpy.linalg.inv(Z))
    da = numpy.zeros(2,float)
    da[0] = numpy.sqrt(da_matrix[0,0])
    da[1] = numpy.sqrt(da_matrix[1,1])

    trueslope = -(beta_k[1]-beta_k[0])

    if (bGraph):
        true = df+trueslope*xaxis 
        fit = a[0] + a[1]*xaxis

        PrintData(xaxis,true,fit,ratio,dratio,'linear')

        name = name + '_' + type + '(linear)'
        PrintPicture(xaxis,true,ratio,dratio,fit,name,figname,'linear')

    if (analytic_uncertainty):
        return a[0],da[0],a[1],da[1]
    else:
        return a[0],a[1]

def SolveNonLin(f,df,a,data,ddata,xaxis,maxiter=20,tol=1e-10):

    K = numpy.size(a)
    usedat = numpy.isfinite(data)
    y = data[usedat]
    nuse = len(y)
    weights = 1.0/ddata[usedat]
    w = numpy.diag(weights) 
    x = xaxis[usedat]
    J = numpy.zeros([nuse,K],dtype = numpy.float64)

    # do the newton-raphson solution
    endnext = False
    for n in range(maxiter):
        
        expt = f(a,x)
        
        J = numpy.transpose(df(a,x))
        WJ = numpy.dot(w,J)
        JTW = numpy.transpose(WJ)
        dy = y - expt
        Z = numpy.dot(JTW,WJ)
        incr_a = numpy.linalg.solve(Z,numpy.dot(JTW,dy))
        a += incr_a
        ra = incr_a/a
        chtol = numpy.sqrt(numpy.dot(ra,ra))
        if (chtol < tol):
            if (endnext == True) or (analytical_estimate == False):
                    endnow == True   # we put in this logic so that we calculate the matrix at the minimum
                                     # if we want the analytic uncertainty 
            endnext = True
            if (endnow == True):
                break

        if (n == maxiter):
             print "Too many iterations for nonlinear least squares"


    da_matrix = numpy.linalg.inv(Z)
    da = numpy.zeros(K,float)
    for k in range(K):
        da[k] = numpy.sqrt(da_matrix[k,k])

    return a,da    

def ExpFit(a,x):
    return numpy.exp(a[0]+a[1]*x)

def dExpFit(a,x):
    e = numpy.exp(a[0]+a[1]*x)
    return [e,x*e]

def NonLinFit(bins,a,N_k,beta_k,U_kn,df=0,analytic_uncertainty=False,bGraph=False,name="",
              figname='nonlin_figure',type='kinetic', tol=1e-10,g=[1,1]):

    # nonlinear model is exp(A + B*E_i), where the i are the bin energies.
    # residuals are y_i - exp(A + B*E_i)
    # dS/dbeta_j = 2\sum_i r_i dr_i/dB_j = 0 
    # 
    # dr_i/dA = exp(A + B*E_i)
    # dr_i/dB = E_i*exp(A + B*E_i)

    K = len(N_k)

    hlist = []
    dhlist = []

    for k in range(0,K):
        hstat = plt.hist(U_kn[k,0:N_k[k]], bins = bins)
        h = (1.0*hstat[0])/N_k[k] 
        hlist.append(h)
        dh = numpy.sqrt(g[k]*h*(1.0-h)/N_k[k])
        dhlist.append(dh)

    ratio = (hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = ratio*(numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2))

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    
    (a,da) = SolveNonLin(ExpFit,dExpFit,a,ratio,dratio,xaxis,tol=tol)

    if (bGraph):
        trueslope = -(beta_k[1]-beta_k[0])
        true = numpy.exp(df+trueslope*xaxis) 
        fit = numpy.exp(a[0] + a[1]*xaxis)

        PrintData(xaxis,true,fit,ratio,dratio,'nonlinear')

        name = name + '_' + type + ' (nonlinear)'
        PrintPicture(xaxis,true,ratio,dratio,fit,name,figname,'nonlinear')

    if (analytic_uncertainty):
        return a[0],da[0],a[1],da[1]
    else:
        return a[0],a[1]

def MaxwellBoltzFit(bins,U,N,kT,figname,name="",type='kinetic',ndof=None,g=1):

    # generate histogram
    hstat = plt.hist(U, bins = bins)
    # normalize the histogram
    h = (1.0*hstat[0])/N 
    # compute the error bars
    dh = numpy.sqrt(g*h*(1.0-h)/N)
    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    

    # we assume we have the correct mean for now, since presumably the temperature works
    mean = numpy.mean(U) 
    std_fit = numpy.std(U)
    std_true = numpy.sqrt(mean*kT)
    if (mean > 25*kT):  #if too big, we use a normal distribution -- we'll use limit of 50 DOF as suggest (by Wikipedia!)
        # note that for a normal distribution, the sample mean and standard deviation give the maximum likelihood information.
        fit = numpy.exp(-(xaxis-mean)**2/(2*std_fit**2))/(numpy.sqrt(2*numpy.pi*std_fit**2))
        true = numpy.exp(-(xaxis-mean)**2/(2*std_true**2))/(numpy.sqrt(2*numpy.pi*std_true**2))
    else:
        # should be a gamma distribution; no std fit
        fit = 2*numpy.sqrt(mean/(numpy.pi*(kT)**3))*exp(-mean/kT)
        if (ndof != None):
            mean_true = 0.5*ndof*kT
            true = 2*numpy.sqrt(meanU/(numpy.pi*(kT)**3))*exp(-meanU/kT)
        else:
            true = None # unless we know the number of DOF, we don't know the true distribution:

    print "--- Kinetic energy analysis ---"
    print ""
    print "kT = %10.4f" % (kT)
    if (ndof == None):
        print "Effective # of DOF = %10.4f" % (2*mean/kT)
    else:     
        print "Reported # of DOF = %10.4f" % ndof 
    if (mean > 25*kT):
        "Approximating the Maxwell-Boltzmann with a normal distribution, as # DOF > 50"
    print "Direct Std = %10.4f, Std from sqrt(U*kT) = %10.4f" % (std_fit,std_true)
    print ""

    name = name + '_' + type + ' str'
    # normalize histogram and error bars
    width = bins[1]-bins[0]  # assumes equal spacing (currently true)
    h /= width
    dh /= width
    PrintPicture(xaxis,true,h,dh,fit,name,figname,'maxwell')
    
def PrintData(xaxis,true,fit,collected,dcollected,type):

    if (type == 'linear'):
        print "----  Linear Fit  ----"
    elif (type == 'nonlinear'):
        print "----  Nonlinear Fit  ----"
    elif (type == 'maxwell'):
        print "----  fit to Maxwell-Boltzmann ----"
    else:
        print "Incorrect type specified!"
        # should die at this point
        return

    print "      X         True     Observed     Error   d(true/obs) sig(true/obs)  Fit   "
    print "---------------------------------------------------------------------------------------"
    for i in range(len(collected)):
        diff = collected[i]-true[i]
        sig = numpy.abs(collected[i]-true[i])/dcollected[i]
        print "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f" % (xaxis[i],true[i],collected[i],dcollected[i],diff,sig,fit[i])


def ProbabilityAnalysis(N_k,T_k,U_kn,kB=0.0083144624,title=None,figname=None,nbins=40,
                        bMaxLikelihood=True,bLinearFit=True,bNonLinearFit=True,reptype=None,nboots=200,g=[1,1],reps=None,cuttails=0.0001,bMaxwell=False):

    K = 2  # should be 2 for now

    beta_k = (1.0/(kB*T_k))

    # determine the bin widths
    maxk = numpy.zeros(K,float)
    mink = numpy.zeros(K,float)

    # cuttails indicates how many we leave out on each tail
    # for now, we choose the range that cuts 0.1% from the tails of the smallest distribution.
    prange = 0.001

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

    if (bMaxwell):
        print "Now fitting to a Maxwell-Boltzmann distribution"        
        for k in range(2):
            fn = figname + '_maxboltz' + str(T_k[k])        
            MaxwellBoltzFit(bins,U_kn[k,0:N_k[k]],N_k[k],kB*T_k[k],fn)

    if (bLinearFit):
        print "Now computing the linear fit parameters"
        fn = figname + '_linear'
        (lindf,dlindf,linslope,dlinslope) = LinFit(bins,N_k,beta_k,U_kn,df=df,name=title,figname=fn,bGraph=True,analytic_uncertainty=True,g=g)
        PrintStats('Linear Fit Analysis (analytical error)',lindf,linslope,kB,T_k,trueslope,ddf=dlindf,dslope=dlinslope)

    if (bNonLinearFit): 
        print "Now computing the nonlinear fit parameters" 
        fn = figname + '_nonlinear'
        (nldf,dnldf,nlslope,dnlslope) = NonLinFit(bins,[df,trueslope],N_k,beta_k,U_kn,df=df,name=title,figname=fn,bGraph=True,analytic_uncertainty=True,g=g)
        PrintStats('Nonlinear Fit Analysis (analytical error)',nldf,nlslope,kB,T_k,trueslope,ddf=dnldf,dslope=dnlslope)

    if (bMaxLikelihood):
        print "Now computing the maximum likelihood parameters" 
        (mldf,dmldf,mlslope,dmlslope) = MaxLikeParams(N_k,beta_k,U_kn,df=df,analytic_uncertainty=True,g=numpy.average(g))
        PrintStats('Maximum Likelihood Analysis (analytical error)',mldf,mlslope,kB,T_k,trueslope,ddf=dmldf,dslope=dmlslope)

    if (reptype == None):
        return

    if (reptype == 'bootstrap'):
        nreps = nboots
        print "Now bootstrapping (n=%d) for uncertainties . . . could take a bit!" % (nboots)
    elif (reptype == 'independent'):
        nreps = len(reps)
        print "Now analyzing %d independent samples . . . could also take a bit!" % (nreps)
    else:
        print "Don't understand reptype = %s; quitting" % (reptype)

    linslopes = numpy.zeros([nreps],float)
    lindfs = numpy.zeros([nreps],float)
    mlslopes = numpy.zeros([nreps],float)
    mldfs = numpy.zeros([nreps],float)
    nlslopes = numpy.zeros([nreps],float)
    nldfs = numpy.zeros([nreps],float)

    if (reptype == 'bootstrap'):    
        Ur_kn = numpy.zeros([K,numpy.max(N_k)],float)

    for n in range(nreps):
        if (n%10 == 0):
            print "Finished %d samples . . ." % (n)
        
        if (reptype == 'bootstrap'):    
            for k in range(K):
                if (g == None):
                    #do normal bootstrap
                    rindex = numpy.random.randint(0,high=N_k[k],size=N_k[k]);  # bootstrap it 
                    Ur_kn[k,0:N_k[k]] = U_kn[k,rindex]
                else:
                    # we are given correlation times.  Do block bootstrapping.
                    gk = int(numpy.ceil(g[k]))
                    nblocks = int(numpy.floor(N_k[k]/gk))
                    # moving blocks bootstrap; all contiguous segments of length gk
                    rindex = numpy.random.randint(0,high=gk*nblocks,size=nblocks); 
                    for nb in range(nblocks):
                        Ur_kn[k,nb*gk:(nb+1)*gk] = U_kn[k,rindex[nb]:rindex[nb]+gk]
                    N_k[k] = nblocks*gk  # we could have a few samples less now
    
        if (reptype == 'independent'):
            Ur_kn = reps[n] 

        if (bLinearFit):    
            (lindfs[n],linslopes[n]) = LinFit(bins,N_k,beta_k,Ur_kn) 

        if (bMaxLikelihood):
            (mldfs[n],mlslopes[n]) = MaxLikeParams(N_k,beta_k,Ur_kn,df=df)

        if (bNonLinearFit):
            (nldfs[n],nlslopes[n]) = NonLinFit(bins,[df,trueslope],N_k,beta_k,Ur_kn,df=df)

    if (bLinearFit):
        PrintStats('Linear Fit Analysis',lindfs,linslopes,kB,T_k,trueslope)

    if (bNonLinearFit):
        PrintStats('Nonlinear Fit Analysis',nldfs,nlslopes,kB,T_k,trueslope)

    if (bMaxLikelihood):
        PrintStats('Maximum Likelihood Analysis',mldfs,mlslopes,kB,T_k,trueslope)

    return
    
