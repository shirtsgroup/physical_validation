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

def check_twodtype(type):  # check if it's a valid type
    if (type=='dbeta-dpressure'):
        print 'Warning: can\'t do 3D fits currently' 
    # throw an exception?
        return False
    else:
        return True

def PrepStrings(type):

    if (type == 'dbeta-constV'):
        pstring = 'ln(P_2(E)/P_1(E))'
        pstringtex = '$\\frac{P_2(E)}{P_1(E)}$'
        pstringlntex = '$\\ln\\frac{P_2(E)}{P_1(E)}$'
        varstring = '$E (kT)$'

    elif (type == 'dbeta-constP'):    
        pstring = 'ln(P_2(H)/P_1(H))'
        pstringtex = '$\\frac{P_2(H)}{P_1(H)}$'
        pstringlntex = '$\\ln\\frac{P_2(H)}{P_1(H)}$'
        varstring = '$H (kT)$'

    elif (type == 'dpressure-constB'):    
        pstring = 'ln(P_2(V)/P_1(V))'
        pstringtex = '$\\frac{P_2(V)}{P_1(V)}$'
        pstringlntex = '$\\ln\\frac{P_2(V)}{P_1(V)}$'
        varstring = '$P (input units)$'

    elif (type == 'dbeta-dpressure'):    
        pstring = ''
        pstringtex = ''
        pstringlntex = ''
        varstring = ''
    else:
        print "Type is not defined!"

    return pstring,pstringtex,pstringlntex,varstring       
        

def PrepInputs(N_k,type='dbeta-constV',beta=None,P=None,U_kn=None,V_kn=None):

    # convenience variables 
    N0 = N_k[0]
    N1 = N_k[1]

    # Currently three types; fitting parameters are: 
    # 1) free energy, dbeta  - constants are beta_ave, variables (vectors) are E 
    # 2) free energy, dpressure - constants are p_ave, variables (vectors) are V  
    # 3) free energy, dbeta, dpressure - constants are beta_ave, p_ave, variables (vectors) are E and V

    if (type == 'dbeta-constV'):
        # allocate space 
        v = numpy.zeros([1,2,numpy.max(N0,N1)],float)        
        vr = numpy.zeros([1,2,numpy.max(N0,N1)],float)
        const = numpy.zeros(1,float)
        dp = numpy.zeros(1,float)

        v[0,0,0:N0] = U_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1])
        dp[0] = beta[0] - beta[1]

    elif (type == 'dbeta-constP'):
        # allocate space 
        v = numpy.zeros([1,2,numpy.max(N0,N1)],float)
        vr = numpy.zeros([1,2,numpy.max(N0,N1)],float)
        const = numpy.zeros(1,float)
        dp[0] = numpy.zeros(1,float)

        v[0,0,0:N0] = U_kn[0,0:N0] + pressure_ave*V_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1] + pressure_ave*V_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1])
        dp = beta[0] - beta[1]
        
    elif (type == 'dpressure'):
        # allocate space 
        v = numpy.zeros([1,2,numpy.max(N0,N1)],float)
        vr = numpy.zeros([1,2,numpy.max(N0,N1)],float)
        const = numpy.zeros(1,float)
        dp = numpy.zeros(1,float)

        v[0,0,0:N0] = V_kn[0,0:N0]
        v[0,1,0:N1] = V_kn[1,0:N1]
        const[0] == 0.5*(P[0] + P[1])
        dp[0] = P[0] - P[1]

    elif (type == 'dbeta-dpressure'):    
        # allocate space 
        v = numpy.zeros([2,2,numpy.max(N0,N1)],float)
        vr = numpy.zeros([2,2,numpy.max(N0,N1)],float)
        const = numpy.zeros(2,float)
        dp = numpy.zeros(2,float)
        v[0,0,0:N0] = U_kn[0,0:N0]
        v[0,1,0:N1] = U_kn[1,0:N1]
        v[1,0,0:N0] = V_kn[0,0:N0]
        v[1,1,0:N1] = V_kn[1,0:N1]
        const[0] = 0.5*(beta[0] + beta[1])
        const[1] = 0.5*(P[0] + P[1])
        dp[0] = beta[0] - beta[1]
        dp[1] = P[0] - P[1]
    else:
        print "Warning:  Type of analysis is not defined!"

    return dp,const,v,vr

#def LogLikelihood(x,N_k,beta_ave,U_kn):
def LogLikelihood(x,N_k,const,v):

    L = len(x)

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log(N1/N0)

    #D0 = M + beta_ave*x[0] + U0*x[1]
    #D1 = M + beta_ave*x[0] + U1*x[1]

    D0 = D1 = M + const[0]*x[0]
    for i in range(L-1):
        D0 = D0 + v[i,0,0:N0]*x[i+1]
        D1 = D1 + v[i,1,0:N1]*x[i+1]

    E0 = 1 + numpy.exp(D0)
    E1 = 1 + numpy.exp(-D1)

    # this is the negative of the log likelihood, since we want to maximize it using fmin

    of = ((numpy.sum(numpy.log(E0)) + numpy.sum(numpy.log(E1))))/N
    #print of

    return of

#def dLogLikelihood(x,N_k,beta,U_kn):
def dLogLikelihood(x,N_k,const,v):
    """
    Derivative with respect to the parameters, to aid the minimization.
    
    """

    L = len(x)

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    M = numpy.log(N1/N0)

    D0 = D1 = M + const[0]*x[0]
    for i in range(L-1):
        D0 = D0 + v[i,0,:]*x[i+1]
        D1 = D1 + v[i,1,:]*x[i+1]

    E0 = 1/(1 + numpy.exp(-D0))
    E1 = 1/(1 + numpy.exp(D1))

    g = numpy.zeros(L,dtype=numpy.float64)

    #this is the gradient of -log likelihood
    #g[0] = (1.0/N)*(numpy.sum(beta*E0) - numpy.sum(beta*E1))
    #g[1] = (1.0/N)*(numpy.sum(U0*E0) - numpy.sum(U1*E1))
    #g[2] = (1.0/N)*(numpy.sum(V0*E0) - numpy.sum(V1*E1))

    g[0] = const[0]*(numpy.sum(E0) - numpy.sum(E1))
    for i in range(L-1):
        g[i+1] = numpy.sum(v[i,0,0:N0]*E0) - numpy.sum(v[i,1,0:N1]*E1)
    return (1.0/N)*g

def d2LogLikelihood(x,N_k,const,v):

    """

    beta = const[0]
    pave = const[1]

    if D = M + beta*x[0] + x[1]*U
    I = \sum_{i=1}^N [[-beta^2/S,-beta*U/S],[-beta*U/S,-U^2/S]] where S = [(1+exp(-D))*(1+exp(D))]

    if D = M + beta*x[0] + x[1]*U + x[2]*V
    I = \sum_{i=1}^N [[-beta^2/S,-beta*U/S,-beta*V/S],[-beta*U/S,-U^2/S,-U*V/S],[-beta*V/S,-U*V^2/S,-V^2/S]] where S = [(1+exp(-D))*(1+exp(D))]
    """

    L = len(x)

    M = numpy.log(N_k[1]/N_k[0])

    N0 = N_k[0]
    N1 = N_k[1]
    N  = N0+N1

    vall = numpy.zeros([L-1,N],dtype=numpy.float64)
    for i in range(L-1):
        vall[i,0:N0] = v[i,0,0:N0]
        vall[i,N0:N] = v[i,1,0:N1]
    
    D = M + const[0]*x[0]
    for i in range(L-1):
        D = D + vall[i,:]*x[i+1]
    
    E = (1 + numpy.exp(-D)) * (1 + numpy.exp(D))

    hf = numpy.zeros([L,L,N],dtype=numpy.float64)

    cones = const[0] * numpy.ones(N,dtype=numpy.float64)

    # fix this to match the     
    for i in range(L):
        if (i == 0):
            a = cones
        else: 
            a = vall[i-1,:]
        for j in range(L):
            if (j == 0): 
                b = cones
            else:
                b = vall[i-1,:]    
            hf[i,j,:] = a*b    

    # this is the hessian of the minimum function (not the max)
    h = -numpy.sum(hf/E,axis=2)/N
                            
    return h

def solveminlike(x, N_k, const,v, tol = 1e-10, maxiter=20):
    
    converge = False
    itol = 1e-2
    rtol = 1e-2
    lasttol = 100000;

    for i in range(maxiter):
        lastx = x
        gx = numpy.transpose(dLogLikelihood(x,N_k,const,v))
        nh = d2LogLikelihood(x,N_k,const,v)
        dx = numpy.linalg.solve(nh,gx)
        x += dx  # adding, not subtracting because of the handling of negatives 
        rx = dx/x
        checktol = numpy.sqrt(numpy.dot(dx,dx))
        checkrtol = numpy.sqrt(numpy.dot(rx,rx))    
        if (checkrtol < tol):
            break
            converge = True
        if (checkrtol > 1.0) and (checktol > lasttol):  # we are possibly diverging. Switch to cg for a bit.
            x = scipy.optimize.fmin_cg(LogLikelihood,lastx,fprime=dLogLikelihood,gtol=itol,args=[N_k,const,v],disp=1)
            itol *= rtol
        lasttol = checktol

    if (i == maxiter) and (converge == False):
        print "Too many iterations, convergence failing"

    return x    

def MaxLikeUncertain(x,N_k,const,v,vave):

    L = len(x)
    d = numpy.zeros(L,float)

    # multiply back by N, since we were dealing with smaller numbers for numerical robustness.
    fi = -(N_k[0] + N_k[1])*d2LogLikelihood(x,N_k,const,v)

    d2 = numpy.linalg.inv(fi)

    # We have a fit to the line y = m(x-Uave) + b, so to get the uncertainty in the free energy back, we need 
    # to add M*Uave back to f.  The uncertainty in cov(b + m*Uave,b+m*Uave) = var(b) + Uave**2*var(m) + Uave*cov(v,m)

    # For two dimensioms, we have the line y = m1(x1-vave1) + m2(x2-vave2) + b
    #  Uncertainty will be cov(b + m1vave1 + m2vave2) = var(b) + vave1^2 var(m1) + vave2^2 var(m2)
    #                                                          + 2vave1 cov(m1,b) + 2vave2 cov(m2,b)
    #                                                          + 2vave1 cov(m1,m2)
    d[0] = const[0]**2*d2[0,0] 
    for i in range(1,L):
        d[0] += vave[i-1]**2*d2[i,i] + 2*vave[i-1]*d2[0,i]   # should this last one be plus or minus
        d[i] = d2[i,i]
        for j in range(i+1,L-1):
            d[0] += 2*vave[i-1]*vave[j-1]*d2[i,j]
    d = numpy.sqrt(d)
    return d

def MaxLikeParams(N_k,dp,const,v,df=0,analytic_uncertainty=False,g=1):

    L = len(const)
    optimum = numpy.zeros(L+1,float)
    vave = numpy.zeros(L,dtype=numpy.float64)
    vmod = numpy.zeros([L,numpy.sum(N_k)],dtype=numpy.float64)
    # for numerical stability, we need to translate the curve
    for i in range(L):
        vave[i] = (numpy.sum(v[i,0,0:N_k[0]]) + numpy.sum(v[i,1,0:N_k[1]]))/numpy.sum(N_k)
        vmod = v - vave

    xstart = numpy.zeros(L+1,float)
    for i in range(L):
        xstart[0] += vave[i]*dp[i]
        xstart[i+1] = dp[i]
    xstart[0] += df
    xstart[0] /= const[0]

    ofit = solveminlike(xstart,N_k,const,vmod,tol=1e-10)

    optimum[0] = ofit[0]*const[0]
    for i in range(L):
        optimum[i+1] = ofit[i+1]
        optimum[0] -= (vave[i]*ofit[i+1])

    if (analytic_uncertainty):
        # incorporate g in an average way.
        doptimum = MaxLikeUncertain(ofit,N_k,const,vmod,vave)*numpy.sqrt(1.0/numpy.average(g))
        return optimum[0], doptimum[0], optimum[1], doptimum[1]
    else:
        return optimum[0], optimum[1]

def Print2DStats(title,type,dfs,slopes,kB,dp,const,trueslope,ddf='N/A',dslope='N/A'):

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
    print "     Estimated slope       vs.   True slope"
    print "---------------------------------------------"
    print "%11.6f +/- %11.6f |   %11.6f" % (slope, dslope, trueslope)
    print "---------------------------------------------"

    quant = numpy.abs((slope-trueslope)/dslope)
    print ""
    print "(That's %.2f quantiles from true slope=%5f, FYI.)" % (quant,trueslope)
    if (quant > 5):
        print " (Ouch!)"
    else:
        print ""

    if (type[0:5] == 'dbeta'):    
        #dp = B1 - B0, const = (B1 + B0)/2, B = 1/kbT
        # so B0 = (const-dp/2), T0 = 1/(kB*B0)
        # so B1 = (const+dp/2), T1 = 1/(kB*B1)
        T0 = (kB*(const-dp/2))**(-1)
        T1 = (kB*(const+dp/2))**(-1)

        print "---------------------------------------------"
        print " True dT = %7.3f, Eff. dT = %7.3f+/-%.3f" % (T0-T1, kB*T0*T1*slope,kB*dslope*T0*T1)
        print "---------------------------------------------"

    if (type == 'dpressure-constT'):
        print "---------------------------------------------"
        print " True dP = %7.3f, Eff. dP = %7.3f+/-%.3f" % (dp, slope, dslope)
        print "---------------------------------------------"

def PrintPicture(xaxis,true,y,dy,fit,type,name,figname,fittype,show=False):

    import matplotlib
    import matplotlib.pyplot as plt

    [pstring,pstringtex,pstringlntex,varstring] = PrepStrings(type)

    print "Now printing figure %s" % (figname)
    plt.clf()
    plt.xlabel(varstring)
    if (fittype == 'linear'):
        plt.title('E vs. log probability ratio for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r"pstringlntex")
        plt.errorbar(xaxis,true,fmt='k-',label = r'$-(\beta_2-\beta_1)E$')
        plt.errorbar(xaxis,fit,fmt='r-',label = 'Fit to $y = b+aB$')
        plt.ylabel(r'$\ln\frac{P_1(E)}{P_2(E)}$')
    elif (fittype == 'nonlinear'):
        plt.title('E vs. probability ratio for ' + name)
        plt.errorbar(xaxis,y,fmt='b-',yerr=dy,label = r"pstringtex")
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

def GenHistogramProbs(N_k,bins,v,g):

    K = len(N_k)
    
    hlist = []
    dhlist = []

    for k in range(0,K):
        hstat = plt.hist(v[0,k,0:N_k[k]], bins = bins)
        h = (1.0*hstat[0])/N_k[k] 
        hlist.append(h)
        dh = numpy.sqrt(g[k]*h*(1.0-h)/N_k[k])
        dhlist.append(dh)

    return hlist,dhlist

def LinFit(bins,N_k,dp,const,v,df=0,analytic_uncertainty=False,bGraph=False,name="",figname='lin_figure',g=[1,1],type='dbeta-constV'):
        
    [hlist,dhlist] = GenHistogramProbs(N_k,bins,v,g)

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

    # the true line is y = df + dp*x, where y is ln P_1(X)/P_2(X)

    if (bGraph):
        trueslope = dp
        true = df+trueslope*xaxis 
        fit = a[0] + a[1]*xaxis

        PrintData(xaxis,true,fit,ratio,dratio,'linear')

        name = name + '(linear)'
        PrintPicture(xaxis,true,ratio,dratio,fit,type,name,figname,'linear')

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

def ExpFit(a,x):  # assume only 2D, since we are not generating histograms
    return numpy.exp(a[0] + a[1]*x)

def dExpFit(a,x):
    s = a[0] + a[1]*x
    e = numpy.exp(s)
    return numpy.array([e,x*e])

def NonLinFit(bins,N_k,dp,const,v,df=0,analytic_uncertainty=False,bGraph=False,name="",
              figname='nonlin_figure', tol=1e-10,g=[1,1], type = 'dbeta-constV'):

    # nonlinear model is exp(A + B*E_i), where the i are the bin energies.
    # residuals are y_i - exp(A + B*E_i)
    # dS/dbeta_j = 2\sum_i r_i dr_i/dB_j = 0 
    # 
    # dr_i/dA = exp(A + B*E_i)
    # dr_i/dB = E_i*exp(A + B*E_i)

    [hlist,dhlist] = GenHistogramProbs(N_k,bins,v,g)

    ratio = (hlist[1]/hlist[0]) # this should have the proper exponential distribution 
    dratio = ratio*(numpy.sqrt((dhlist[0]/hlist[0])**2 + (dhlist[1]/hlist[1])**2))

    xaxis = (bins[0:len(bins)-1] + bins[1:len(bins)])/2    

    # starting point for nonlinear fit
    a = numpy.zeros(len(dp)+1)
    a[0] = df
    a[1:len(dp)+1] = dp[:]

    (a,da) = SolveNonLin(ExpFit,dExpFit,a,ratio,dratio,xaxis,tol=tol)

    if (bGraph):
        trueslope = dp
        true = numpy.exp(df+trueslope*xaxis) 
        fit = ExpFit(a,xaxis)
        
        PrintData(xaxis,true,fit,ratio,dratio,'nonlinear')

        name = name + ' (nonlinear)'
        PrintPicture(xaxis,true,ratio,dratio,fit,type,name,figname,'nonlinear')

    if (analytic_uncertainty):
        return a[0],da[0],a[1],da[1]
    else:
        return a[0],a[1]

def MaxwellBoltzFit(bins,U,N,kT,figname,name="",ndof=None,g=1):

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

    # name = name + str
    # normalize histogram and error bars
    width = bins[1]-bins[0]  # assumes equal spacing (currently true)
    h /= width
    dh /= width
    PrintPicture(xaxis,true,h,dh,fit,'dbeta-constV',name,figname,'maxwell')
    
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


def ProbabilityAnalysis(N_k,type='dbeta-constV',T_k=None,P_k=None,U_kn=None,V_kn=None,kB=0.0083144624,title=None,figname=None,nbins=40,
                        bMaxLikelihood=True,bLinearFit=True,bNonLinearFit=True,reptype=None,nboots=200,g=[1,1],reps=None,cuttails=0.0001,bMaxwell=False):

    K = len(N_k)  # should be 2 pretty much always

    # initialize constant terms
    beta_ave = None
    P_ave = None

    if (T_k != None):
        beta_k = (1.0/(kB*T_k))
        beta_ave = numpy.average(beta_k)
    if (P_k != None):    
        P_ave = numpy.average(P_k)

    # turn the prepare the variables we are going to work with    
    [dp,const,v,vr] = PrepInputs(N_k,type,beta_k,P_k,U_kn,V_kn)
    [pstring,pstringtex,pstringlntex,varstring] = PrepStrings(type)
    
    if (check_twodtype(type)):  # if it's 2D, we can graph, otherwise, there is too much histogram error 
        # determine the bin widths
        maxk = numpy.zeros(K,float)
        mink = numpy.zeros(K,float)

        # cuttails indicates how many we leave out on each tail
        # for now, we choose the range that cuts 0.1% from the tails of the smallest distribution.
        prange = 0.001

        for k in range(K):
            maxk[k] = scipy.stats.scoreatpercentile(v[0,k,0:N_k[k]],100*(1-prange))
            mink[k] = scipy.stats.scoreatpercentile(v[0,k,0:N_k[k]],100*(prange))

        binmax = numpy.min(maxk)
        binmin = numpy.max(mink)

        bins = numpy.zeros(nbins+1,float)
        for i in range(nbins+1):
            bins[i] = binmin + (binmax-binmin)*(i/(1.0*nbins))    

    #===================================================================================================
    # Calculate free energies with different methods
    #===================================================================================================    

    if ((type == 'dbeta-constV') or (type == 'dbeta-constP')):         
        trueslope = (beta_k[0]-beta_k[1])
        w_F = (beta_k[1]-beta_k[0])*v[0,0,0:N_k[0]];
        w_R = (beta_k[0]-beta_k[1])*v[0,1,0:N_k[1]];

        print "True slope of %s should be %.8f" % (pstring,trueslope)

    if (type == 'dpressure-constT'):         
        trueslope = (P_k[0]-P_k[1])

        w_F = (P_k[1]-P_k[0])*beta_ave*v[0,0,0:N_k[0]];
        w_R = (P_k[0]-P_k[1])*beta_ave*v[0,1,0:N_k[1]];

        print "True slope of %s should be %.8f" % (pstring,trueslope)

    if (type == 'dbeta-dpressure'):
        w_F =  (beta_k[1]-beta_k[0])*v[0,0,0:N_k[0]] + (beta_k[1]*p_K[1]-beta_k[0]*P_k[0])*v[1,0,0:N_k[0]]
        w_R =  (beta_k[0]-beta_k[1])*v[0,1,0:N_k[1]] + (beta_k[0]*p_K[0]-beta_k[1]*P_k[1])*v[1,1,0:N_k[1]]

    print "Now computing log of partition functions using BAR"
    
    (df,ddf) = pymbar.BAR(w_F,w_R,compute_uncertainty=True)

    print "using %.5f for log of partition functions computed from BAR" % (df) 
    print "Uncertainty in quantity is %.5f" % (ddf)
    print "Assuming this is negligible compared to sampling error at individual points" 

    if (bMaxwell):  # only applies for kinetic energies
        print "Now fitting to a Maxwell-Boltzmann distribution"        
        for k in range(2):
            fn = figname + '_maxboltz' + str(T_k[k])        
            MaxwellBoltzFit(bins,U_kn[k,0:N_k[k]],N_k[k],kB*T_k[k],fn)

    if (bLinearFit and check_twodtype(type)):
        print "Now computing the linear fit parameters"
        fn = figname + '_linear'
        (lindf,dlindf,linslope,dlinslope) = LinFit(bins,N_k,dp,const,v,df=df,name=title,figname=fn,bGraph=True,analytic_uncertainty=True,g=g)
        Print2DStats('Linear Fit Analysis (analytical error)',type,lindf,linslope,kB,dp,const,trueslope,ddf=dlindf,dslope=dlinslope)

    if (bNonLinearFit and check_twodtype(type)): 
        print "Now computing the nonlinear fit parameters" 
        fn = figname + '_nonlinear'
        (nldf,dnldf,nlslope,dnlslope) = NonLinFit(bins,N_k,dp,const,v,df=df,name=title,figname=fn,bGraph=True,analytic_uncertainty=True,g=g)
        Print2DStats('Nonlinear Fit Analysis (analytical error)',type,nldf,nlslope,kB,dp,const,trueslope,ddf=dnldf,dslope=dnlslope)

    if (bMaxLikelihood):
        print "Now computing the maximum likelihood parameters" 
        (mldf,dmldf,mlslope,dmlslope) = MaxLikeParams(N_k,dp,const,v,df=df,analytic_uncertainty=True,g=numpy.average(g))
        if (check_twodtype(type)):
            Print2DStats('Maximum Likelihood Analysis (analytical error)',type,mldf,mlslope,kB,dp,const,trueslope,ddf=dmldf,dslope=dmlslope)

    if (reptype == None):
        return

    if (reptype == 'bootstrap'):
        nreps = nboots
        print "Now bootstrapping (n=%d) for uncertainties . . . could take a bit of time!" % (nboots)
    elif (reptype == 'independent'):
        nreps = len(reps)
        print "Now analyzing %d independent samples . . . could take a bit of time!" % (nreps)
    else:
        print "Don't understand reptype = %s; quitting" % (reptype)

    linslopes = numpy.zeros([nreps],float)
    lindfs = numpy.zeros([nreps],float)
    mlslopes = numpy.zeros([nreps],float)
    mldfs = numpy.zeros([nreps],float)
    nlslopes = numpy.zeros([nreps],float)
    nldfs = numpy.zeros([nreps],float)

    for n in range(nreps):
        if (n%10 == 0):
            print "Finished %d samples . . ." % (n)
        
        if (reptype == 'bootstrap'):    
            for k in range(K):
                if (g == None):
                    #do normal bootstrap
                    rindex = numpy.random.randint(0,high=N_k[k],size=N_k[k]);  # bootstrap it 
                    for i in range(len(const)):
                        vr[i,k,0:N_k[k]] = v[i,k,rindex]
                else:
                    # we are given correlation times.  Do block bootstrapping.
                    gk = int(numpy.ceil(g[k]))
                    nblocks = int(numpy.floor(N_k[k]/gk))
                    # moving blocks bootstrap; all contiguous segments of length gk
                    rindex = numpy.random.randint(0,high=gk*nblocks,size=nblocks); 
                    for nb in range(nblocks):
                        for i in range(len(const)):
                            vr[i,k,nb*gk:(nb+1)*gk] = v[i,k,rindex[nb]:rindex[nb]+gk]
                    N_k[k] = nblocks*gk  # we could have a few samples less now
    
        if (reptype == 'independent'):
            vr = reps[n] 

        if (bLinearFit):    
            (lindfs[n],linslopes[n]) = LinFit(bins,N_k,dp,const,vr) 

        if (bNonLinearFit):
            (nldfs[n],nlslopes[n]) = NonLinFit(bins,N_k,dp,const,vr,df=df)

        if (bMaxLikelihood):
            (mldfs[n],mlslopes[n]) = MaxLikeParams(N_k,dp,const,vr,df=df)

    if (check_twodtype):    
        if (bLinearFit):
            Print2DStats('Linear Fit Analysis',type,lindfs,linslopes,kB,dp,const,trueslope)

        if (bNonLinearFit):
            Print2DStats('Nonlinear Fit Analysis',type,nldfs,nlslopes,kB,dp,const,trueslope)

        if (bMaxLikelihood):
            Print2DStats('Maximum Likelihood Analysis',type,mldfs,mlslopes,kB,dp,const,trueslope)

    return
    
