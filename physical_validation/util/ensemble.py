###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
This file largely corresponds to the checkensemble.py code originally
published on https://github.com/shirtsgroup/checkensemble. It now serves
as the low-level functionality of the high-level module
:mod:`physical_validation.ensemble`.
"""
from __future__ import division
import numpy as np
import scipy.optimize

import pymbar

from . import trajectory
from . import error as pv_error
from . import plot


def generate_histograms(traj1, traj2, g1, g2, bins):

    n1 = np.size(traj1)
    n2 = np.size(traj2)

    h1, _ = np.histogram(traj1, bins=bins, density=True)
    h2, _ = np.histogram(traj2, bins=bins, density=True)
    dh1 = np.sqrt(g1 * h1 * (1 - h1) / n1)
    dh2 = np.sqrt(g2 * h2 * (1 - h2) / n2)

    # used to be:
    h1, _ = np.histogram(traj1, bins=bins)
    h1 = h1/n1
    dh1 = np.sqrt(g1 * h1 * (1 - h1) / n1)
    h2, _ = np.histogram(traj2, bins=bins)
    h2 = h2/n2
    dh2 = np.sqrt(g2 * h2 * (1 - h2) / n2)

    return h1, h2, dh1, dh2


def do_linear_fit(traj1, traj2, g1, g2, bins,
                  screen=False, filename=None,
                  trueslope=0.0, trueoffset=0.0,
                  name=None, units=None):

    h1, h2, dh1, dh2 = generate_histograms(traj1, traj2, g1, g2, bins)

    #  v  copied from checkensemble.py  v
    ratio = np.log(h2 / h1)
    dratio = np.sqrt((dh1/h1)**2 + (dh2/h2)**2)

    usedat = np.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0/dratio[usedat]

    xaxis = (bins[:-1] + bins[1:])/2
    x = xaxis[usedat]

    x_mat = np.ones([nuse, 2])
    x_mat[:, 1] = x

    w = np.diag(weights) 
    wx = np.dot(w, x_mat)
    wy = np.dot(w, y)
    wx_t = np.transpose(wx)
    z = np.dot(wx_t, wx)
    wxy = np.dot(wx_t, wy)

    a = np.linalg.solve(z, wxy)
    da_matrix = np.transpose(np.linalg.inv(z))
    da = np.zeros(2)
    da[0] = np.sqrt(da_matrix[0, 0])
    da[1] = np.sqrt(da_matrix[1, 1])

    # the true line is y = df + dp*x, where y is ln P_1(X)/P_2(X)
    #  ^  end copied from checkensemble.py  ^

    do_plot = screen or filename is not None
    if do_plot:
        true = trueoffset+trueslope*xaxis
        fit = a[0] + a[1]*xaxis

        data = [{'x': xaxis,
                 'y': ratio,
                 'y_err': dratio,
                 'name': 'Simulation'},
                {'x': xaxis,
                 'y': fit,
                 'name': 'Fit to simulation'},
                {'x': xaxis,
                 'y': true,
                 'name': 'Analytical ratio'}]

        if units is not None:
            units = ' [' + units + ']'
        else:
            units = ''

        annot = ('{:.1f}'.format(abs((a[1] - trueslope) / da[1])) +
                 ' quantiles')

        plot.plot(data,
                  legend='best',
                  title='Log probability ratio',
                  xlabel='Energy' + units,
                  ylabel=r'$\log\frac{P_2(E)}{P_1(E)}$',
                  filename=filename,
                  screen=screen,
                  axtext=annot)

    return a, da


def do_max_likelihood_fit(traj1, traj2, g1, g2,
                          init_params=None):

    # ============================================================= #
    # Define (negative) log-likelihood function and its derivatives #
    # ============================================================= #
    def log_likelihood(a, ene1, ene2):
        # Returns negative of eq (8) of check_ensemble paper
        #
        # Uses log (1/f(x)) == -log(f(x))
        # and log(1 + e^x) == log(e^x (e^-x + 1)) == x + log(1 + e^-x)
        #     ^(a)                                   ^(b)
        # form (a) -> 0 for x->-inf, -> inf for x->inf
        # form (b) -> NaN for x->-inf, -> x for x->inf
        # combined: -> 0 for x-> -inf, -> x for x-> inf
        def log_1_plus_exp(y):
            def f(yy):
                with np.errstate(over='raise'):
                    try:
                        xx = np.log(1 + np.exp(yy))
                    except FloatingPointError:
                        xx = yy + np.log(1 + np.exp(-yy))
                    return xx
            return np.vectorize(f)(y)

        if a.size == 2:
            return (np.sum(log_1_plus_exp(a[0] + a[1]*ene1)) +
                    np.sum(log_1_plus_exp(-a[0] - a[1]*ene2)))
        else:
            return (np.sum(log_1_plus_exp(a[0] + a[1]*ene1[0] + a[2]*ene1[1])) +
                    np.sum(log_1_plus_exp(-a[0] - a[1]*ene2[0] - a[2]*ene2[1])))

    def da_log_likelihood(a, ene1, ene2):
        # Returns the first derivative wrt the parameters a of log_likelihood
        #
        # d/da0 log(1 + exp(a0 + a1*E)) == exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == 1 / (1 + exp(-a0 - a1*E))
        # d/da1 log(1 + exp(a0 + a1*E)) == E * exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == E / (1 + exp(-a0 - a1*E))
        def inv_1_plus_exp(y):
            def f(yy):
                with np.errstate(over='raise'):
                    try:
                        xx = 1. / (1 + np.exp(yy))
                    except FloatingPointError:
                        xx = 0.
                    return xx
            return np.vectorize(f)(y)

        if a.size == 2:
            d = np.zeros(2)
            d[0] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1)) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2)))
            d[1] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1) * ene1) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2) * ene2))
        else:
            d = np.zeros(3)
            d[0] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1])) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1])))
            d[1] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1]) * ene1[0]) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1]) * ene2[0]))
            d[2] = (np.sum(inv_1_plus_exp(-a[0] - a[1]*ene1[0] - a[2]*ene1[1]) * ene1[1]) -
                    np.sum(inv_1_plus_exp(a[0] + a[1]*ene2[0] + a[2]*ene2[1]) * ene2[1]))

        return d

    def hess_log_likelihood(a, ene1, ene2):
        # Returns the hessian wrt the parameters a of log_likelihood
        # fac1 = 1 / (2 + 2*cosh(a0 + a1*ene1))
        # h1 = [[ fac1,      ene1*fac1    ],
        #       [ ene1*fac1, ene1**2*fac1 ]]
        # fac2 = 1 / (2 + 2*cosh(a0 + a1*ene2))
        # h2 = [[ fac2,      ene2*fac2    ],
        #       [ ene2*fac2, ene2**2*fac2 ]]
        # h = h1 + h2

        if a.size == 2:
            fac1 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene1))
            fac2 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene2))

            h = np.zeros((2, 2))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[0, 1] = h[1, 0] = np.sum(ene1 * fac1) + np.sum(ene2 * fac2)
            h[1, 1] = np.sum(ene1 * ene1 * fac1) + np.sum(ene2 * ene2 * fac2)

        else:
            fac1 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene1[0] + a[2]*ene1[1]))
            fac2 = 1 / (2 + 2*np.cosh(a[0] + a[1]*ene2[0] + a[2]*ene2[1]))

            h = np.zeros((3, 3))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[1, 1] = np.sum(ene1[0] * ene1[0] * fac1) + np.sum(ene2[0] * ene2[0] * fac2)
            h[2, 2] = np.sum(ene1[1] * ene1[1] * fac1) + np.sum(ene2[1] * ene2[1] * fac2)

            h[0, 1] = h[1, 0] = np.sum(ene1[0] * fac1) + np.sum(ene2[0] * fac2)
            h[0, 2] = h[2, 0] = np.sum(ene1[1] * fac1) + np.sum(ene2[1] * fac2)
            h[1, 2] = h[2, 1] = (np.sum(ene1[0] * ene1[1] * fac1) +
                                 np.sum(ene2[0] * ene2[1] * fac2))

        return h

    # ==================================================== #
    # Minimize the negative of the log likelihood function #
    # ==================================================== #
    if init_params is None:
        init_params = np.zeros(traj1.ndim + 1)

    min_res = scipy.optimize.minimize(
        log_likelihood,
        x0=init_params,
        args=(traj1, traj2),
        method='slsqp',
        jac=da_log_likelihood
    )

    # fallback options
    if not min_res.success:
        prev_method = 'slsqp'
        for method in ['nelder-mead']:
            print('Note: Max-Likelihood minimization failed using {:s} method. '
                  'Trying {:s}.'.format(prev_method, method))
            min_res = scipy.optimize.minimize(
                log_likelihood,
                x0=init_params,
                args=(traj1, traj2),
                method=method
            )
            if min_res.success:
                break
            else:
                prev_method = method
        if not min_res.success:
            for method in ['bfgs', 'cg']:
                print('Note: Max-Likelihood minimization failed using {:s} method. '
                      'Trying {:s}.'.format(prev_method, method))
                min_res = scipy.optimize.minimize(
                    log_likelihood,
                    x0=init_params,
                    args=(traj1, traj2),
                    method=method,
                    jac=da_log_likelihood
                )
                if min_res.success:
                    break
                else:
                    prev_method = method

        if not min_res.success:
            for method in ['dogleg', 'trust-ncg']:
                print('Note: Max-Likelihood minimization failed using {:s} method. '
                      'Trying {:s}.'.format(prev_method, method))
                min_res = scipy.optimize.minimize(
                    log_likelihood,
                    x0=init_params,
                    args=(traj1, traj2),
                    method=method,
                    jac=da_log_likelihood,
                    hess=hess_log_likelihood
                )
                if min_res.success:
                    break

    final_params = min_res.x

    # ======================= #
    # Calculate uncertainties #
    # ======================= #
    cov = np.linalg.inv(hess_log_likelihood(final_params, traj1, traj2))
    final_error = np.sqrt(np.diag(cov))*np.sqrt(np.average([g1, g2]))

    return final_params, final_error


def check_bins(traj1, traj2, bins):
    # check for empty bins
    h1, _ = np.histogram(traj1, bins=bins)
    h2, _ = np.histogram(traj2, bins=bins)
    empty = np.where((h1 == 0) | (h2 == 0))[0]

    if np.size(empty) == 0:
        return bins
    elif np.size(empty) == 1:
        empty = empty[0]
        if empty > np.size(bins) / 2:
            return bins[:empty]
        else:
            return bins[empty+1:]
    else:
        # find longest non-empty interval
        max_interval = np.argmax(empty[1:] - empty[:-1])
        left = empty[max_interval] + 1
        right = empty[max_interval + 1]
        return bins[left:right]


def print_stats(title,
                fitvals, dfitvals,
                kb, param1, param2, trueslope,
                dtemp=False, dpress=False, dmu=False,
                dtempdpress=False, dtempdmu=False):

    # if simple 1d:
    #     fitvals = [df, slope]
    #     dfitvals = [ddf, dslope]
    # if simple 2d:
    #     fitvals = [df, slope0, slope1]
    #     dfitvals = [ddf, dslope0, dslope1]
    # if bootstrapped 1d:
    #     fitvals = [[df, slope], [df, slope], ...]
    #     dfitvals = None
    # if bootstrapped 2d:
    #     fitvals = [[df, slope0, slope1], [df, slope0, slope1], ...]
    #     dfitvals = None
    if fitvals.ndim > 1:
        dfitvals = np.std(fitvals, axis=0)
        fitvals = np.average(fitvals, axis=0)

    if np.ndim(trueslope) == 0:
        trueslopes = np.array([trueslope])
    else:
        trueslopes = trueslope

    free_energy = fitvals[0]
    slopes = fitvals[1:]
    dfree_energy = dfitvals[0]
    dslopes = dfitvals[1:]

    print('='*50)
    print(title)
    print('='*50)
    print('Free energy')
    print('    {:.5f} +/- {:.5f}'.format(free_energy, dfree_energy))
    print('{:27s}      |  {:s}'.format('Estimated slope', 'True slope'))
    for slope, dslope, trueslope in zip(slopes, dslopes, trueslopes):
        print('    {:<9.6f} +/- {:<9.6f}      |  {:<9.6f}'.format(slope, dslope, trueslope))
        quant = np.abs((slope-trueslope)/dslope)
        print('    ({:.2f} quantiles from true slope)'.format(quant))

    if dtemp or dtempdpress or dtempdmu:
        # slope is estimated beta2 - beta1
        # kb * slope == 1/T1' - 1/T2' == (T2' - T1')/(T1'*T2')
        # So we'll assume dT' == T1' - T2' ~= kb * slope * T1*T2
        slope = slopes[0]
        dslope = dslopes[0]
        if dtemp:
            t1 = param1
            t2 = param2
        else:
            t1 = param1[0]
            t2 = param2[0]
        print('{:27s}      |  {:s}'.format('Estimated dT', 'True dT'))
        print('    {:<6.3f} +/- {:<6.3f}            |  {:<6.3f}'.format(
            kb * slope * t1 * t2,
            kb * dslope * t1 * t2,
            t1 - t2
        ))
    if dpress or dtempdpress:
        # slope is estimated P2 - P1
        if dpress:
            slope = slopes[0]
            dslope = dslopes[0]
            trueslope = trueslopes[0]
        else:
            slope = slopes[1]
            dslope = dslopes[1]
            trueslope = trueslopes[1]
        print('{:27s}      |  {:s}'.format('Estimated dP', 'True dP'))
        print('    {:<6.3f} +/- {:<6.3f}            |  {:<6.3f}'.format(
            -slope, np.abs(dslope), -trueslope
        ))
    print('='*50)


def check_1d(traj1, traj2, param1, param2, kb,
             quantity, dtemp=False, dpress=False, dmu=False,
             nbins=40, cutoff=0.001, seed=None,
             verbosity=1, screen=False, filename=None):
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1 : array-like
        Trajectory of the first simulation
        If dtemp:
            NVT: Potential energy U or total energy E = U + K
            NPT: Enthalpy H = U + pV or total energy E = H + K
        If dpress:
            NPT: Volume V
    traj2 : array-like
        Trajectory of the second simulation
        If dtemp:
            NVT: Potential energy U or total energy E = U + K
            NPT: Enthalpy H = U + pV or total energy E = H + K
        If dpress:
            NPT: Volume V
    param1 : float
        Target temperature or pressure of the first simulation
    param2 : float
        Target temperature or pressure of the second simulation
    kb : float
        Boltzmann constant in same units as the energy trajectories
    nbins : int
        Number of bins used to assess distributions of the trajectories
        Default: 40
    cutoff : float
        Tail cutoff of distributions.
        Default: 0.001 (0.1%)
    seed : int
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping non-reproducible.
    verbosity : int
        Verbosity level.
        Default: 1 (only most important output)
    screen : bool, optional
        Plot distributions on screen.
        Default: False.
    filename : string, optional
        Plot distributions to `filename`.pdf.
        Default: None.

    Returns
    -------

    """

    if (not (dtemp or dpress or dmu) or
       (dtemp and dpress) or
       (dtemp and dmu) or
       (dpress and dmu)):
        raise pv_error.InputError(['dtemp', 'dpress', 'dmu'],
                                  'Need to specify exactly one of `dtemp`, `dpress` and `dmu`.')

    if dmu:
        raise NotImplementedError('check_1d: Testing of `dmu` not implemented.')

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = 'ln(P_2(' + quantity + ')/P_1(' + quantity + '))'
    trueslope = 0
    if dtemp:
        trueslope = 1/(kb * param1) - 1/(kb * param2)
    elif dpress:
        trueslope = param2 - param1

    if verbosity > 1:
        print('Analytical slope of {:s}: {:.8f}'.format(pstring, trueslope))

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.equilibrate(traj1, verbose=(verbosity > 1), name='Trajectory 1')
    traj1 = trajectory.decorrelate(traj1, verbose=(verbosity > 1), name='Trajectory 1')
    traj2 = trajectory.equilibrate(traj2, verbose=(verbosity > 1), name='Trajectory 2')
    traj2 = trajectory.decorrelate(traj2, verbose=(verbosity > 1), name='Trajectory 2')

    # calculate inefficiency
    g1 = pymbar.timeseries.statisticalInefficiency(traj1)
    g2 = pymbar.timeseries.statisticalInefficiency(traj2)

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full, traj2=traj2_full, cut=cutoff
    )
    if not min_ene:
        raise pv_error.InputError(['traj1', 'traj2'],
                                  'No overlap between trajectories.')
    if verbosity > 0:
        print('Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.'.format(
            traj1.shape[0] / traj1_full.shape[0],
            traj2.shape[0] / traj2_full.shape[0]
        ))

    # calculate bins
    bins = np.linspace(min_ene, max_ene, nbins+1)
    bins = check_bins(traj1, traj2, bins)
    if np.size(bins) < 3:
        raise pv_error.InputError(['traj1', 'traj2', 'nbins', 'cutoff'],
                                  'Less than 3 bins were filled in the overlap region.\n'
                                  'Ensure sufficient overlap between the trajectories, and '
                                  'consider increasing `cutoff` or `nbins` if there is '
                                  'sufficient overlap but unusually long tails.')

    w_f = -trueslope * traj1_full
    w_r = trueslope * traj2_full

    if verbosity > 2:
        print('Computing log of partition functions using pymbar.BAR...')
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print('Using {:.5f} for log of partition functions as computed from BAR.'.format(df))
        print('Uncertainty in quantity is {:.5f}.'.format(ddf))
        print('Assuming this is negligible compared to sampling error at individual points.')

    # ========== #
    # linear fit #
    # ========== #
    if verbosity > 2:
        print('Computing linear fit parameters (for plotting / comparison)')

    fitvals, dfitvals = do_linear_fit(
        traj1=traj1_full, traj2=traj2_full, g1=g1, g2=g2, bins=bins,
        screen=screen, filename=filename,
        trueslope=trueslope, trueoffset=df,
        name=None, units=None
    )

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant['linear'] = [abs((slope - trueslope)/dslope)]
    if verbosity > 1:
        print_stats(
            title='Linear Fit Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            dtemp=dtemp, dpress=dpress, dmu=dmu
        )

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print('Computing the maximum likelihood parameters')

    fitvals, dfitvals = do_max_likelihood_fit(traj1_full, traj2_full, g1, g2)

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant['maxLikelihood'] = [abs((slope - trueslope)/dslope)]
    if verbosity > 0:
        print_stats(
            title='Maximum Likelihood Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            dtemp=dtemp, dpress=dpress, dmu=dmu
        )

    return quant['maxLikelihood']


def check_2d(traj1, traj2, param1, param2, kb, pvconvert,
             quantity, dtempdpress=False, dtempdmu=False,
             nbins=40, cutoff=0.001, seed=None,
             verbosity=1, screen=False, filename=None):
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1 : array-like, 2d
        Trajectory of the first simulation
        If dtempdpress:
            traj[0,:]: Potential energy U or total energy E = U + K
            traj[1,:]: Volume V
    traj2 : array-like, 2d
        Trajectory of the second simulation
        If dtempdpress:
            traj[0,:]: Potential energy U or total energy E = U + K
            traj[1,:]: Volume V
    param1 : array-like
        If dtempdpress:
            Target temperature and pressure of the first simulation
    param2 : array-like
        If dtempdpress:
            Target temperature and pressure of the first simulation
    kb : float
        Boltzmann constant in same units as the energy trajectories
    pvconvert : float
        Conversion from pressure * volume to energy units
    nbins : int
        Number of bins used to assess distributions of the trajectories
        Default: 40
    cutoff : float
        Tail cutoff of distributions.
        Default: 0.001 (0.1%)
    seed : int
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping non-reproducible.
    verbosity : int
        Verbosity level.
        Default: 1 (only most important output)
    screen : bool, optional
        Plot distributions on screen.
        Default: False.
    filename : string, optional
        Plot distributions to `filename`.pdf.
        Default: None.

    Returns
    -------

    """

    if not (dtempdpress or dtempdmu) or (dtempdpress and dtempdmu):
        raise pv_error.InputError(['dtempdpress', 'dtempdmu'],
                                  'Need to specify exactly one of `dtempdpress` and `dtempdmu`.')

    if dtempdmu:
        raise NotImplementedError('check_2d: Testing of `dtempdmu` not implemented.')

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = ('ln(P_2(' + quantity[0] + ', ' + quantity[1] + ')/' +
               'P_1(' + quantity[0] + ', ' + quantity[1] + '))')
    trueslope = np.zeros(2)
    facs = [None, None]
    if dtempdpress:
        trueslope = np.array([
            1/(kb * param1[0]) - 1/(kb * param2[0]),
            pvconvert*(1/(kb * param1[0]) * param1[1] - 1/(kb * param2[0]) * param2[1])
        ])
        facs = [[1, param1[1]], [1, param2[1]]]

    if verbosity > 1:
        print('Analytical slope of {:s}: {:.8f}, {:.8f}'.format(
            pstring, trueslope[0], trueslope[1]
        ))

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.equilibrate(traj1, verbose=(verbosity > 1), name='Trajectory 1')
    traj1 = trajectory.decorrelate(traj1, facs=facs[0], verbose=(verbosity > 1), name='Trajectory 1')
    traj2 = trajectory.equilibrate(traj2, verbose=(verbosity > 1), name='Trajectory 2')
    traj2 = trajectory.decorrelate(traj2, facs=facs[1], verbose=(verbosity > 1), name='Trajectory 2')

    # calculate inefficiency
    g1 = np.array([
            pymbar.timeseries.statisticalInefficiency(traj1[0]),
            pymbar.timeseries.statisticalInefficiency(traj1[1])
        ])
    g2 = np.array([
            pymbar.timeseries.statisticalInefficiency(traj2[0]),
            pymbar.timeseries.statisticalInefficiency(traj2[1])
        ])

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full, traj2=traj2_full, cut=cutoff
    )
    if min_ene is None:
        raise pv_error.InputError(['traj1', 'traj2'],
                                  'No overlap between trajectories.')
    if verbosity > 0:
        print('Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.'.format(
            traj1.shape[1] / traj1_full.shape[1],
            traj2.shape[1] / traj2_full.shape[1]
        ))

    w_f = -trueslope[0] * traj1_full[0] - trueslope[1] * traj1_full[1]
    w_r = trueslope[0] * traj2_full[0] + trueslope[1] * traj2_full[1]

    if verbosity > 2:
        print('Computing log of partition functions using pymbar.BAR...')
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print('Using {:.5f} for log of partition functions as computed from BAR.'.format(df))
        print('Uncertainty in quantity is {:.5f}.'.format(ddf))
        print('Assuming this is negligible compared to sampling error at individual points.')

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print('Computing the maximum likelihood parameters')

    fitvals, dfitvals = do_max_likelihood_fit(traj1_full, traj2_full, g1, g2)

    slope = fitvals[1:]
    dslope = dfitvals[1:]
    quant['maxLikelihood'] = np.abs((slope - trueslope)/dslope)
    if verbosity > 0:
        print_stats(
            title='Maximum Likelihood Analysis (analytical error)',
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            dtempdpress=dtempdpress, dtempdmu=dtempdmu
        )

    return quant['maxLikelihood']
