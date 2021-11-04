###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Pascal T. Merz <pascal.merz@me.com>                       #
#               Michael R. Shirts <michael.shirts@colorado.edu>           #
#                                                                         #
#    Copyright (c) 2017-2021 University of Colorado Boulder               #
#              (c) 2012      The University of Virginia                   #
#                                                                         #
###########################################################################
r"""
This file reimplements most functionality of the checkensemble.py code
originally published on https://github.com/shirtsgroup/checkensemble. It
serves as the low-level functionality of the high-level module
:mod:`physical_validation.ensemble`.
"""
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pymbar
import scipy.optimize

from . import error as pv_error
from . import plot, trajectory


def chemical_potential_energy(
    chemical_potential: np.ndarray,
    number_of_species: np.ndarray,
) -> np.ndarray:
    r"""
    Calculates the chemical potential energy of a trajectory by
    returning the sum of the current number of species multiplied
    by the chemical potential.
    """
    assert chemical_potential.ndim == 1
    assert number_of_species.ndim == 1 or number_of_species.ndim == 2
    if number_of_species.ndim == 1:
        assert chemical_potential.size == 1
        return chemical_potential[0] * number_of_species
    else:
        assert number_of_species.shape[1] == chemical_potential.size
        # chemical_potential: 1 x num_pot
        # number_of_species:  num_frames x num_pot
        return np.sum(chemical_potential * number_of_species, axis=1)


def generate_histograms(
    traj1: np.ndarray, traj2: np.ndarray, g1: float, g2: float, bins: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    n1 = np.size(traj1)
    n2 = np.size(traj2)

    h1 = np.histogram(traj1, bins=bins)[0] / n1
    h2 = np.histogram(traj2, bins=bins)[0] / n2
    dh1 = np.sqrt(g1 * h1 * (1 - h1) / n1)
    dh2 = np.sqrt(g2 * h2 * (1 - h2) / n2)

    return h1, h2, dh1, dh2


def do_linear_fit(
    traj1: np.ndarray,
    traj2: np.ndarray,
    g1: float,
    g2: float,
    bins: np.ndarray,
    screen: bool,
    filename: Optional[str],
    trueslope: float,
    trueoffset: float,
    units: Optional[str],
    xlabel: str,
    ylabel: str,
) -> Tuple[np.ndarray, np.ndarray]:

    h1, h2, dh1, dh2 = generate_histograms(traj1, traj2, g1, g2, bins)

    #  v  copied from checkensemble.py  v
    ratio = np.log(h2 / h1)
    dratio = np.sqrt((dh1 / h1) ** 2 + (dh2 / h2) ** 2)

    usedat = np.isfinite(ratio)
    y = ratio[usedat]
    nuse = len(y)
    weights = 1.0 / dratio[usedat]

    xaxis = (bins[:-1] + bins[1:]) / 2
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
        true = trueoffset + trueslope * xaxis
        fit = a[0] + a[1] * xaxis

        data = [
            {"x": xaxis, "y": ratio, "y_err": dratio, "name": "Simulation"},
            {"x": xaxis, "y": fit, "name": "Fit to simulation"},
            {"x": xaxis, "y": true, "name": "Analytical ratio"},
        ]

        if units is not None:
            units = " [" + units + "]"
        else:
            units = ""

        annot = (
            "{:.1f}".format(abs((a[1] - trueslope) / da[1]))
            + " quantiles (linear estimate)"
        )

        plot.plot(
            data,
            legend="upper left",
            title="Log probability ratio",
            xlabel=xlabel + units,
            ylabel=ylabel,
            filename=filename,
            screen=screen,
            axtext=annot,
        )

    return a, da


def do_max_likelihood_fit(
    traj1: np.ndarray,
    traj2: np.ndarray,
    g1: Union[float, np.ndarray],
    g2: Union[float, np.ndarray],
    init_params: np.ndarray,
    verbose: bool,
) -> Tuple[np.ndarray, np.ndarray]:

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
                with np.errstate(over="raise"):
                    try:
                        xx = np.log(1 + np.exp(yy))
                    except FloatingPointError:
                        xx = yy + np.log(1 + np.exp(-yy))
                    return xx

            return np.vectorize(f)(y)

        if a.size == 2:
            return np.sum(log_1_plus_exp(a[0] + a[1] * ene1)) + np.sum(
                log_1_plus_exp(-a[0] - a[1] * ene2)
            )
        else:
            return np.sum(
                log_1_plus_exp(a[0] + a[1] * ene1[0] + a[2] * ene1[1])
            ) + np.sum(log_1_plus_exp(-a[0] - a[1] * ene2[0] - a[2] * ene2[1]))

    def da_log_likelihood(a, ene1, ene2):
        # Returns the first derivative wrt the parameters a of log_likelihood
        #
        # d/da0 log(1 + exp(a0 + a1*E)) == exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == 1 / (1 + exp(-a0 - a1*E))
        # d/da1 log(1 + exp(a0 + a1*E)) == E * exp(a0 + a1*E) / (1 + exp(a0 + a1*E))
        #                               == E / (1 + exp(-a0 - a1*E))
        def inv_1_plus_exp(y):
            def f(yy):
                with np.errstate(over="raise"):
                    try:
                        xx = 1.0 / (1 + np.exp(yy))
                    except FloatingPointError:
                        xx = 0.0
                    return xx

            return np.vectorize(f)(y)

        if a.size == 2:
            d = np.zeros(2)
            d[0] = np.sum(inv_1_plus_exp(-a[0] - a[1] * ene1)) - np.sum(
                inv_1_plus_exp(a[0] + a[1] * ene2)
            )
            d[1] = np.sum(inv_1_plus_exp(-a[0] - a[1] * ene1) * ene1) - np.sum(
                inv_1_plus_exp(a[0] + a[1] * ene2) * ene2
            )
        else:
            d = np.zeros(3)
            d[0] = np.sum(
                inv_1_plus_exp(-a[0] - a[1] * ene1[0] - a[2] * ene1[1])
            ) - np.sum(inv_1_plus_exp(a[0] + a[1] * ene2[0] + a[2] * ene2[1]))
            d[1] = np.sum(
                inv_1_plus_exp(-a[0] - a[1] * ene1[0] - a[2] * ene1[1]) * ene1[0]
            ) - np.sum(inv_1_plus_exp(a[0] + a[1] * ene2[0] + a[2] * ene2[1]) * ene2[0])
            d[2] = np.sum(
                inv_1_plus_exp(-a[0] - a[1] * ene1[0] - a[2] * ene1[1]) * ene1[1]
            ) - np.sum(inv_1_plus_exp(a[0] + a[1] * ene2[0] + a[2] * ene2[1]) * ene2[1])

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
            fac1 = 1 / (2 + 2 * np.cosh(a[0] + a[1] * ene1))
            fac2 = 1 / (2 + 2 * np.cosh(a[0] + a[1] * ene2))

            h = np.zeros((2, 2))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[0, 1] = h[1, 0] = np.sum(ene1 * fac1) + np.sum(ene2 * fac2)
            h[1, 1] = np.sum(ene1 * ene1 * fac1) + np.sum(ene2 * ene2 * fac2)

        else:
            fac1 = 1 / (2 + 2 * np.cosh(a[0] + a[1] * ene1[0] + a[2] * ene1[1]))
            fac2 = 1 / (2 + 2 * np.cosh(a[0] + a[1] * ene2[0] + a[2] * ene2[1]))

            h = np.zeros((3, 3))

            h[0, 0] = np.sum(fac1) + np.sum(fac2)
            h[1, 1] = np.sum(ene1[0] * ene1[0] * fac1) + np.sum(
                ene2[0] * ene2[0] * fac2
            )
            h[2, 2] = np.sum(ene1[1] * ene1[1] * fac1) + np.sum(
                ene2[1] * ene2[1] * fac2
            )

            h[0, 1] = h[1, 0] = np.sum(ene1[0] * fac1) + np.sum(ene2[0] * fac2)
            h[0, 2] = h[2, 0] = np.sum(ene1[1] * fac1) + np.sum(ene2[1] * fac2)
            h[1, 2] = h[2, 1] = np.sum(ene1[0] * ene1[1] * fac1) + np.sum(
                ene2[0] * ene2[1] * fac2
            )

        return h

    # ==================================================== #
    # Minimize the negative of the log likelihood function #
    # ==================================================== #
    min_res = checkensemble_solver(
        fun=log_likelihood,
        x0=init_params,
        args=(traj1, traj2),
        jac=da_log_likelihood,
        hess=hess_log_likelihood,
    )

    # fallback options
    if not min_res.success:
        # built-in was unsuccessful
        if verbose:
            print(
                "Note: Max-Likelihood minimization failed using built-in method. "
                "Trying scipy method 'dogleg'."
            )
        for variation in [1.0, 0.9, 1.1]:
            min_res = scipy.optimize.minimize(
                log_likelihood,
                x0=init_params * variation,
                args=(traj1, traj2),
                method="dogleg",
                jac=da_log_likelihood,
                hess=hess_log_likelihood,
                options=dict(initial_trust_radius=0.001),
            )
            if min_res.success:
                break

    if not min_res.success:
        # dogleg was unsuccessful
        if verbose:
            print(
                "Note: Max-Likelihood minimization failed using built-in method. "
                "Trying scipy method 'nelder-mead'."
            )
        for variation in [1.0, 0.9, 1.1]:
            min_res = scipy.optimize.minimize(
                log_likelihood,
                x0=init_params * variation,
                args=(traj1, traj2),
                method="nelder-mead",
            )
            if min_res.success:
                break

    if not min_res.success:
        raise RuntimeError("MaxLikelihood: Unable to minimize function.")

    final_params = min_res.x

    # ======================= #
    # Calculate uncertainties #
    # ======================= #
    cov = np.linalg.inv(hess_log_likelihood(final_params, traj1, traj2))
    final_error = np.sqrt(np.diag(cov)) * np.sqrt(np.average([g1, g2]))

    return final_params, final_error


def checkensemble_solver(
    fun, x0, args, jac, hess, tol=1e-10, maxiter=20
) -> scipy.optimize.OptimizeResult:
    # This is the solver used in checkensemble, unchanged except for some
    # modernization / adaptation to coding style
    success = False
    x = np.array(x0)
    itol = 1e-2
    rtol = 1e-2
    lasttol = float("inf")
    lastx = x

    it = 0
    gx = 0
    nh = 0
    for it in range(maxiter):
        gx = np.transpose(jac(x, *args))
        nh = hess(x, *args)
        dx = np.linalg.solve(nh, gx)
        x -= dx
        rx = dx / x
        checktol = np.sqrt(np.dot(dx, dx))
        checkrtol = np.sqrt(np.dot(rx, rx))

        if checkrtol < tol:
            success = True
            break

        if checkrtol > 1.0 and checktol > lasttol:
            x = scipy.optimize.fmin_cg(
                fun, lastx, fprime=jac, gtol=itol, args=args, disp=False
            )
            itol *= rtol

        lasttol = checktol
        lastx = x

    return scipy.optimize.OptimizeResult(
        x=x,
        success=success,
        nit=it,
        status=int(not success),
        fun=fun(x, *args),
        jac=gx,
        hess=nh,
        nfev=int(np.round(np.log(itol * 100) / np.log(rtol))),
        njev=it,
        nhev=it,
    )


def check_bins(traj1: np.ndarray, traj2: np.ndarray, bins: np.ndarray) -> np.ndarray:
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
            return bins[empty + 1 :]
    else:
        # find longest non-empty interval
        empty = np.insert(np.append(empty, [40]), 0, [-1])
        max_interval = np.argmax(empty[1:] - empty[:-1])
        left = empty[max_interval] + 1
        right = empty[max_interval + 1]
        return bins[left:right]


def print_stats(
    title: str,
    fitvals: np.ndarray,
    dfitvals: Optional[np.ndarray],
    kb: float,
    param1: Union[float, np.ndarray],
    param2: Union[float, np.ndarray],
    trueslope: Union[float, np.ndarray],
    temp: Optional[float],
    pvconvert: Optional[float],
    dtemp: bool,
    dpress: bool,
    dmu: bool,
    dtempdpress: bool,
    dtempdmu: bool,
) -> None:
    # if simple 1d:
    #     fitvals = [df, slope]
    #     dfitvals = [ddf, dslope]
    # if simple 2d:
    #     fitvals = [df, slope0, slope1]
    #     dfitvals = [ddf, dslope0, dslope1]
    # if bootstrapped 1d:
    #     fitvals = [[df, slope], [df_i, slope_i], ...]
    #     dfitvals = None
    # if bootstrapped 2d:
    #     fitvals = [[df, slope0, slope1], [df_i, slope0_i, slope1_i], ...]
    #     dfitvals = None

    if fitvals.ndim > 1:  # bootstrapped
        slopes = fitvals[0, 1:]
        dslopes = np.std(fitvals[1:, 1:], axis=0)
        free_energy = fitvals[0, 0]
        dfree_energy = np.std(fitvals[1:, 0], axis=0)
    else:
        slopes = fitvals[1:]
        free_energy = fitvals[0]
        if dfitvals is not None:
            dslopes = dfitvals[1:]
            dfree_energy = dfitvals[0]
        else:
            dslopes = np.zeros(slopes.size)
            dfree_energy = np.zeros(1)

    if np.ndim(trueslope) == 0:
        trueslopes = np.array([trueslope])
    else:
        trueslopes = trueslope

    print("=" * 50)
    print(title)
    print("=" * 50)
    print("Free energy")
    print("    {:.5f} +/- {:.5f}".format(free_energy, dfree_energy))
    print("{:27s}      |  {:s}".format("Estimated slope", "True slope"))
    for slope, dslope, trueslope in zip(slopes, dslopes, trueslopes):
        print(
            "    {:<9.6f} +/- {:<9.6f}      |  {:<9.6f}".format(
                slope, dslope, trueslope
            )
        )
        quant = np.abs((slope - trueslope) / dslope)
        print("    ({:.2f} quantiles from true slope)".format(quant))

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
        print("{:27s}      |  {:s}".format("Estimated dT", "True dT"))
        print(
            "    {:<6.1f} +/- {:<6.1f}            |  {:<6.1f}".format(
                kb * slope * t1 * t2, kb * dslope * t1 * t2, t2 - t1
            )
        )
    if dpress or dtempdpress:
        # slope is estimated (P1 - P2)/beta*pvconvert (1d), or
        #                    (P1/b1 - P2/b2)*pvconvert (2d)
        # in 2d, we'll assume dP == P1 - P2 ~= slope * .5(T1 + T2) / pvconvert
        if temp is None and dtempdpress:
            temp = 0.5 * (param1[0] + param2[0])
        if dpress:
            press = -slopes[0] * (kb * temp) / pvconvert
            ddpress = -dslopes[0] * (kb * temp) / pvconvert
            truepress = -trueslopes[0] * (kb * temp) / pvconvert
        else:
            press = -slopes[1] * (kb * temp) / pvconvert
            ddpress = -dslopes[1] * (kb * temp) / pvconvert
            truepress = -trueslopes[1] * (kb * temp) / pvconvert
        print("{:27s}      |  {:s}".format("Estimated dP", "True estimated dP"))
        print(
            "    {:<6.1f} +/- {:<6.1f}            |  {:<6.1f}".format(
                press, np.abs(ddpress), truepress
            )
        )
    if dmu or dtempdmu:
        # slope is estimated (mu1 - mu2)/beta (1d), or
        #                    (mu1/b1 - mu2/b2) (2d)
        # in 2d, we'll assume dmu == mu1 - mu2 ~= slope * .5(mu1 + mu2)
        if temp is None and dtempdmu:
            temp = 0.5 * (param1[0] + param2[0])
        if dmu:
            idx = 0
        else:
            idx = 1

        mu = -slopes[idx] * (kb * temp)
        ddmu = -dslopes[idx] * (kb * temp)
        truemu = -trueslopes[idx] * (kb * temp)
        print("{:27s}      |  {:s}".format("Estimated dmu", "True estimated dmu"))
        print(
            "    {:<6.1f} +/- {:<6.1f}            |  {:<6.1f}".format(
                mu, np.abs(ddmu), truemu
            )
        )
    print("=" * 50)


def estimate_interval(
    ens_string: str,
    ens_temp: float,
    energy: np.ndarray,
    kb: float,
    ens_press: Optional[float],
    volume: Optional[np.ndarray],
    pvconvert: Optional[float],
    ens_mu: Optional[np.ndarray],
    species_number: Optional[np.ndarray],
    verbosity: int,
    cutoff: float,
    tunit: str,
    punit: str,
    munit: str,
    data_is_uncorrelated: bool,
) -> Dict[str, Union[float, List[float]]]:
    result = {}
    if ens_string == "NVT":
        # Discard burn-in period and time-correlated frames
        energy = trajectory.prepare(
            energy,
            cut=cutoff,
            verbosity=verbosity - 1,
            name="Energy",
            skip_preparation=data_is_uncorrelated,
        )
        # dT
        sig = np.std(energy)
        result["dT"] = 2 * kb * ens_temp * ens_temp / sig
    elif ens_string == "NPT":
        enthalpy = energy + pvconvert * ens_press * volume
        traj_2d = np.array([energy, volume])
        # Discard burn-in period and time-correlated frames
        enthalpy = trajectory.prepare(
            enthalpy,
            cut=cutoff,
            verbosity=verbosity - 1,
            name="Enthalpy",
            skip_preparation=data_is_uncorrelated,
        )
        volume_1d = trajectory.prepare(
            volume,
            cut=cutoff,
            verbosity=verbosity - 1,
            name="Volume",
            skip_preparation=data_is_uncorrelated,
        )
        traj_2d = trajectory.prepare(
            traj_2d,
            cut=cutoff,
            verbosity=verbosity - 1,
            name="2D-Trajectory",
            skip_preparation=data_is_uncorrelated,
        )

        # dT
        sig = np.std(enthalpy)
        result["dT"] = 2 * kb * ens_temp * ens_temp / sig
        # dP
        sig = np.std(volume_1d) * pvconvert
        result["dP"] = 2 * kb * ens_temp / sig
        # dTdP
        cov = np.cov(traj_2d)
        sig = np.sqrt(np.diag(cov))
        sig[1] *= pvconvert
        result["dTdP"] = [
            2 * kb * ens_temp * ens_temp / sig[0],
            2 * kb * ens_temp / sig[1],
        ]
    elif ens_string == "muVT":
        if ens_mu.size > 1:
            print(
                "Note: Your muVT ensemble has more than one mu value. Ensemble check is "
                "implemented for\n"
                "      * ensembles differing in T only, with any number of mu values\n"
                "      * ensembles differing in one value of mu only\n"
                "      * ensembles differing in T and one value of mu.\n"
            )
        else:
            species_number = species_number.flatten()
        muvt_energy = energy - chemical_potential_energy(ens_mu, species_number)
        # Discard burn-in period and time-correlated frames
        muvt_energy = trajectory.prepare(
            muvt_energy,
            cut=cutoff,
            verbosity=verbosity - 1,
            name="muVT Energy",
            skip_preparation=data_is_uncorrelated,
        )
        try:
            species_number_traj = trajectory.prepare(
                species_number,
                cut=cutoff,
                verbosity=verbosity - 1,
                name="N",
                skip_preparation=data_is_uncorrelated,
            )
        except NotImplementedError:
            raise NotImplementedError(
                "Preparing of the trajectory is not implemented for more than two dimensions. "
                "Your trajectory has 3 or more species. If you are sure that your trajectory "
                "is uncorrelated, you can retry this command using the argument "
                "`data_is_uncorrelated=True`."
            )

        traj_2d = np.array([energy, species_number]) if ens_mu.size == 1 else None
        if traj_2d is not None:
            traj_2d = trajectory.prepare(
                traj_2d,
                cut=cutoff,
                verbosity=verbosity - 1,
                name="2D-Trajectory",
                skip_preparation=data_is_uncorrelated,
            )

        # dT
        sig = np.std(muvt_energy)
        result["dT"] = 2 * kb * ens_temp * ens_temp / sig
        # dmu
        sig = np.array([np.std(species_number_traj, axis=0)]).flatten()
        if sig.size == 1:
            result["dmu"] = 2 * kb * ens_temp / sig[0]
        else:
            result["dmu"] = [2 * kb * ens_temp / s for s in sig]
        # dTdP
        if traj_2d is not None:
            cov = np.cov(traj_2d)
            sig = np.sqrt(np.diag(cov))
            result["dTdmu"] = [
                2 * kb * ens_temp * ens_temp / sig[0],
                2 * kb * ens_temp / sig[1],
            ]
    else:
        raise pv_error.InputError("ens_str", "Unrecognized ensemble string.")

    if verbosity > 0:
        print(
            "A rule of thumb states that good error recognition can be expected when\n"
            "spacing the tip of the distributions by about two standard deviations.\n"
            "Based on this rule, and the assumption that the standard deviation of the\n"
            "distributions is largely independent of the state point, here's an estimate\n"
            "for the interval given the current simulation:"
        )
        if ens_string == "NVT":
            print("Current trajectory: NVT, T = {:.2f} {:s}".format(ens_temp, tunit))
            print("Suggested interval: dT = {:.1f} {:s}".format(result["dT"], tunit))
        if ens_string == "NPT":
            print(
                "Current trajectory: NPT, T = {:.2f} {:s}, P = {:.2f} {:s}".format(
                    ens_temp, tunit, ens_press, punit
                )
            )
            print("Suggested interval:")
            print("  Temperature-only: dT = {:.1f} {:s}".format(result["dT"], tunit))
            print("  Pressure-only: dP = {:.1f} {:s}".format(result["dP"], punit))
            print(
                "  Combined: dT = {:.1f} {:s}, dP = {:.1f} {:s}".format(
                    result["dTdP"][0], tunit, result["dTdP"][1], punit
                )
            )
        if ens_string == "muVT":

            def print_mu(mu, format_string):
                mu = np.array([mu]).flatten()
                if mu.size > 1:
                    return np.array2string(
                        mu, formatter={"float_kind": lambda x: format_string.format(x)}
                    )
                else:
                    return format_string.format(mu[0])

            print(
                "Current trajectory: muVT, T = {:.2f} {:s}, mu = {:s} {:s}".format(
                    ens_temp, tunit, print_mu(ens_mu, "{:.2f}"), munit
                )
            )
            print("Suggested interval:")
            print("  Temperature-only: dT = {:.1f} {:s}".format(result["dT"], tunit))
            print(
                "  Chemical potential-only: dmu = {:s} {:s}".format(
                    print_mu(result["dmu"], "{:.1f}"), munit
                )
            )
            if "dTdmu" in result:
                print(
                    "  Combined: dT = {:.1f} {:s}, dmu = {:.1f} {:s}".format(
                        result["dTdmu"][0], tunit, result["dTdmu"][1], munit
                    )
                )

    return result


def check_1d(
    traj1: np.ndarray,
    traj2: np.ndarray,
    param1: float,
    param2: float,
    kb: float,
    quantity: str,
    dtemp: bool,
    dpress: bool,
    dmu: bool,
    temp: Optional[float],
    pvconvert: Optional[float],
    nbins: int,
    cutoff: float,
    bootstrap_seed: Optional[int],
    bootstrap_error: bool,
    bootstrap_repetitions: int,
    verbosity: int,
    screen: bool,
    filename: Union[str],
    xlabel: str,
    xunit: Optional[str],
    data_is_uncorrelated: bool,
) -> List[float]:
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1
        Trajectory of the first simulation
        If dtemp:

            * NVT: Potential energy U or total energy E = U + K
            * NPT: Enthalpy H = U + pV or total energy E = H + K

        If dpress:

            * NPT: Volume V

    traj2
        Trajectory of the second simulation
        If dtemp:

            * NVT: Potential energy U or total energy E = U + K
            * NPT: Enthalpy H = U + pV or total energy E = H + K

        If dpress:

            * NPT: Volume V

    param1
        Target temperature or pressure of the first simulation
    param2
        Target temperature or pressure of the second simulation
    kb
        Boltzmann constant in same units as the energy trajectories
    quantity
        Name of quantity analyzed (used for printing only)
    dtemp
        Set to True if trajectories were simulated at different temperature
    dpress
        Set to True if trajectories were simulated at different pressure
    temp
        The temperature in equal temperature, differing pressure NPT simulations.
        Needed to print optimal dP.
    pvconvert
        Conversion from pressure * volume to energy units.
        Needed to print optimal dP.
    dmu
        Set to True if trajectories were simulated at different chemical potential
    nbins
        Number of bins used to assess distributions of the trajectories
    cutoff
        Tail cutoff of distributions.
    bootstrap_seed
        Sets the random number seed for bootstrapping.
        If set, bootstrapping will be reproducible.
        If `None`, bootstrapping is non-reproducible.
    bootstrap_error
        Calculate the standard error via bootstrap resampling
    bootstrap_repetitions
        Number of bootstrap repetitions drawn
    verbosity
        Verbosity level.
    screen
        Plot distributions on screen.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    xlabel
        x-axis label used for plotting
    xunit
        x-axis label unit used for plotting
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
        The number of quantiles the computed result is off the analytical one.

    """

    if (
        not (dtemp or dpress or dmu)
        or (dtemp and dpress)
        or (dtemp and dmu)
        or (dpress and dmu)
    ):
        raise pv_error.InputError(
            ["dtemp", "dpress", "dmu"],
            "Need to specify exactly one of `dtemp`, `dpress` and `dmu`.",
        )

    if dpress and (temp is None or pvconvert is None):
        raise pv_error.InputError(
            ["dpress", "temp", "pvconvert"],
            "`ensemble.check_1d` with `dpress=True` requires `temp` and `pvconvert`.",
        )

    if dmu and temp is None:
        raise pv_error.InputError(
            ["dmu", "temp"], "`ensemble.check_1d` with `dmu=True` requires `temp`."
        )

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = "ln(P_2(" + quantity + ")/P_1(" + quantity + "))"
    trueslope = 0
    if dtemp:
        trueslope = 1 / (kb * param1) - 1 / (kb * param2)
    elif dpress:
        trueslope = (param1 - param2) / (kb * temp) * pvconvert
    elif dmu:
        trueslope = -(param1 - param2) / (kb * temp)

    if verbosity > 1:
        print("Analytical slope of {:s}: {:.8f}".format(pstring, trueslope))

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.prepare(
        traj1,
        cut=cutoff,
        verbosity=verbosity,
        name="Trajectory 1",
        skip_preparation=data_is_uncorrelated,
    )
    traj2 = trajectory.prepare(
        traj2,
        cut=cutoff,
        verbosity=verbosity,
        name="Trajectory 2",
        skip_preparation=data_is_uncorrelated,
    )

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full, traj2=traj2_full
    )
    if verbosity > 0:
        print(
            "Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.".format(
                traj1.shape[0] / traj1_full.shape[0],
                traj2.shape[0] / traj2_full.shape[0],
            )
        )
    if verbosity > 0 and dtemp:
        sig1 = np.std(traj1_full)
        sig2 = np.std(traj2_full)
        dt1 = 2 * kb * param1 * param1 / sig1
        dt2 = 2 * kb * param2 * param2 / sig2
        if verbosity > 1:
            print(
                "A rule of thumb states that a good overlap is found when dT/T = (2*kB*T)/(sig),\n"
                "where sig is the standard deviation of the energy distribution.\n"
                "For the current trajectories, dT = {:.1f}, sig1 = {:.1f} and sig2 = {:.1f}.\n"
                "According to the rule of thumb, given T1, a good dT is dT = {:.1f}, and\n"
                "                                given T2, a good dT is dT = {:.1f}.".format(
                    param2 - param1, sig1, sig2, dt1, dt2
                )
            )
        print(
            "Rule of thumb estimates that dT = {:.1f} would be optimal "
            "(currently, dT = {:.1f})".format(0.5 * (dt1 + dt2), param2 - param1)
        )
    if verbosity > 0 and dpress:
        sig1 = np.std(traj1_full) * pvconvert
        sig2 = np.std(traj2_full) * pvconvert
        dp1 = 2 * kb * temp / sig1
        dp2 = 2 * kb * temp / sig2
        if verbosity > 1:
            print(
                "A rule of thumb states that a good overlap is found when dP = (2*kB*T)/(sig),\n"
                "where sig is the standard deviation of the volume distribution.\n"
                "For the current trajectories, dP = {:.1f}, sig1 = {:.1g} and sig2 = {:.1g}.\n"
                "According to the rule of thumb, given P1, a good dP is dP = {:.1f}, and\n"
                "                                given P2, a good dP is dP = {:.1f}.".format(
                    param2 - param1, sig1, sig2, dp1, dp2
                )
            )
        print(
            "Rule of thumb estimates that dP = {:.1f} would be optimal "
            "(currently, dP = {:.1f})".format(0.5 * (dp1 + dp2), param2 - param1)
        )
    if verbosity > 0 and dmu:
        sig1 = np.std(traj1_full)
        sig2 = np.std(traj2_full)
        dmu1 = 2 * kb * temp / sig1
        dmu2 = 2 * kb * temp / sig2
        if verbosity > 1:
            print(
                "A rule of thumb states that a good overlap is found when dP = (2*kB*T)/(sig),\n"
                "where sig is the standard deviation of the volume distribution.\n"
                "For the current trajectories, dmu = {:.1f}, sig1 = {:.1g} and sig2 = {:.1g}.\n"
                "According to the rule of thumb, given mu1, a good dmu is dmu = {:.1f}, and\n"
                "                                given mu2, a good dmu is dmu = {:.1f}.".format(
                    param2 - param1, sig1, sig2, dmu1, dmu2
                )
            )
        print(
            "Rule of thumb estimates that mu = {:.1f} would be optimal "
            "(currently, dmu = {:.1f})".format(
                0.5 * (dmu1 + dmu2), abs(param2 - param1)
            )
        )
    if not min_ene:
        raise pv_error.InputError(
            ["traj1", "traj2"], "No overlap between trajectories."
        )
    # calculate bins
    bins = np.linspace(min_ene, max_ene, nbins + 1)
    bins = check_bins(traj1, traj2, bins)
    if np.size(bins) < 3:
        raise pv_error.InputError(
            ["traj1", "traj2", "nbins", "cutoff"],
            "Less than 3 bins were filled in the overlap region.\n"
            "Ensure sufficient overlap between the trajectories, and "
            "consider increasing `cutoff` or `nbins` if there is "
            "sufficient overlap but unusually long tails.",
        )

    # calculate inefficiency
    g1 = pymbar.timeseries.statisticalInefficiency(traj1)
    g2 = pymbar.timeseries.statisticalInefficiency(traj2)

    w_f = -trueslope * traj1
    w_r = trueslope * traj2

    if verbosity > 2:
        print("Computing log of partition functions using pymbar.BAR...")
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print(
            "Using {:.5f} for log of partition functions as computed from BAR.".format(
                df
            )
        )
        print("Uncertainty in quantity is {:.5f}.".format(ddf))
        print(
            "Assuming this is negligible compared to sampling error at individual points."
        )

    # ========== #
    # linear fit #
    # ========== #
    if verbosity > 2:
        print("Computing linear fit parameters (for plotting / comparison)")

    fitvals, dfitvals = do_linear_fit(
        traj1=traj1,
        traj2=traj2,
        g1=g1,
        g2=g2,
        bins=bins,
        screen=screen,
        filename=filename,
        trueslope=trueslope,
        trueoffset=df,
        units=xunit,
        xlabel=xlabel,
        ylabel=r"$\log\frac{P_2(" + quantity + ")}{P_1(" + quantity + ")}$",
    )

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant["linear"] = [abs((slope - trueslope) / dslope)]
    if verbosity > 1:
        print_stats(
            title="Linear Fit Analysis (analytical error)",
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=temp,
            pvconvert=pvconvert,
            dtemp=dtemp,
            dpress=dpress,
            dmu=dmu,
            dtempdmu=False,
            dtempdpress=False,
        )

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print("Computing the maximum likelihood parameters")

    fitvals, dfitvals = do_max_likelihood_fit(
        traj1,
        traj2,
        g1,
        g2,
        init_params=np.array([df, trueslope]),
        verbose=(verbosity > 1),
    )

    slope = fitvals[1]
    dslope = dfitvals[1]
    quant["maxLikelihood"] = [abs((slope - trueslope) / dslope)]
    if (verbosity > 0 and not bootstrap_error) or verbosity > 1:
        print_stats(
            title="Maximum Likelihood Analysis (analytical error)",
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=temp,
            pvconvert=pvconvert,
            dtemp=dtemp,
            dpress=dpress,
            dmu=dmu,
            dtempdmu=False,
            dtempdpress=False,
        )

    if not bootstrap_error:
        return quant["maxLikelihood"]

    # =============================== #
    # bootstrapped max-likelihood fit #
    # =============================== #
    if verbosity > 0:
        print(
            "Computing bootstrapped maximum likelihood parameters... "
            "[0/{:d}]".format(bootstrap_repetitions),
            end="",
        )

    if bootstrap_seed is not None:
        np.random.seed(bootstrap_seed)

    bs_fitvals = []
    for n, (t1, t2) in enumerate(
        zip(
            trajectory.bootstrap(traj1, bootstrap_repetitions),
            trajectory.bootstrap(traj2, bootstrap_repetitions),
        )
    ):

        # use overlap region
        t1, t2, min_ene, max_ene = trajectory.overlap(traj1=t1, traj2=t2)
        # calculate inefficiency
        g1 = pymbar.timeseries.statisticalInefficiency(t1)
        g2 = pymbar.timeseries.statisticalInefficiency(t2)
        # calculate max_likelihood fit
        fv, _ = do_max_likelihood_fit(
            t1,
            t2,
            g1,
            g2,
            init_params=np.array([df, trueslope]),
            verbose=(verbosity > 2),
        )
        bs_fitvals.append(fv)
        # print progress
        if verbosity > 0:
            print(
                "\rComputing bootstrapped maximum likelihood parameters... "
                "[{:d}/{:d}]".format(n + 1, bootstrap_repetitions),
                end="",
            )

    print()
    bs_fitvals = np.array(bs_fitvals)
    # slope = np.average(fitvals[:, 1])
    dslope = np.std(bs_fitvals[:, 1], axis=0)
    quant["bootstrap"] = [abs((slope - trueslope) / dslope)]
    if verbosity > 0:
        print_stats(
            title="Maximum Likelihood Analysis (bootstrapped error)",
            fitvals=np.concatenate(([fitvals], bs_fitvals)),
            dfitvals=None,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=temp,
            pvconvert=pvconvert,
            dtemp=dtemp,
            dpress=dpress,
            dmu=dmu,
            dtempdmu=False,
            dtempdpress=False,
        )

    return quant["bootstrap"]


def check_2d(
    traj1: np.ndarray,
    traj2: np.ndarray,
    param1: np.ndarray,
    param2: np.ndarray,
    kb: float,
    pvconvert: Optional[float],
    quantity: List[str],
    dtempdpress: bool,
    dtempdmu: bool,
    cutoff: float,
    bootstrap_seed: Optional[int],
    bootstrap_error: bool,
    bootstrap_repetitions: int,
    verbosity: int,
    screen: bool,
    filename: Optional[str],
    data_is_uncorrelated: bool,
) -> List[float]:
    r"""
    Checks whether the energy trajectories of two simulation performed at
    different temperatures have sampled distributions at the analytically
    expected ratio.

    Parameters
    ----------
    traj1
        Trajectory of the first simulation
        If dtempdpress:

            * traj[0,:]: Potential energy U or total energy E = U + K
            * traj[1,:]: Volume V
    traj2
        Trajectory of the second simulation
        If dtempdpress:

            * traj[0,:]: Potential energy U or total energy E = U + K
            * traj[1,:]: Volume V
    param1
        If dtempdpress:
            Target temperature and pressure of the first simulation
    param2
        If dtempdpress:
            Target temperature and pressure of the first simulation
    kb
        Boltzmann constant in same units as the energy trajectories
    pvconvert
        Conversion from pressure * volume to energy units
    quantity
        Names of quantities analyzed (used for printing only)
    dtempdpress
        Set to True if trajectories were simulated at different
        temperature and pressure
    dtempdmu
        Set to True if trajectories were simulated at different
        temperature and chemical potential
    cutoff
        Tail cutoff of distributions.
    bootstrap_seed
        Sets the random number seed for bootstrapping.
        If set, bootstrapping will be reproducible.
        If `None`, bootstrapping is non-reproducible.
    bootstrap_error
        Calculate the standard error via bootstrap resampling
    bootstrap_repetitions
        Number of bootstrap repetitions drawn
    verbosity
        Verbosity level.
    screen
        Plot distributions on screen.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.
        Default: False

    Returns
    -------
        The number of quantiles the computed result is off the analytical one.

    """

    if not (dtempdpress or dtempdmu) or (dtempdpress and dtempdmu):
        raise pv_error.InputError(
            ["dtempdpress", "dtempdmu"],
            "Need to specify exactly one of `dtempdpress` and `dtempdmu`.",
        )

    if screen or filename is not None:
        raise NotImplementedError("check_2d: Plotting not implemented.")

    if dtempdpress and pvconvert is None:
        raise pv_error.InputError(
            "pvconvert", "When using `dtempdpress`, `pvconvert` is required."
        )

    # =============================== #
    # prepare constants, strings etc. #
    # =============================== #
    pstring = (
        "ln(P_2("
        + quantity[0]
        + ", "
        + quantity[1]
        + ")/"
        + "P_1("
        + quantity[0]
        + ", "
        + quantity[1]
        + "))"
    )

    factor = -1
    if dtempdpress:
        factor = pvconvert

    # The true slope in the first dimension is always beta_1 - beta_2,
    # in the second dimension it's either pvconvert * (beta_1*P_1 - beta_2*P_2)
    # or -(beta_1*mu_1 - beta_2*mu_2)
    trueslope = np.array(
        [
            1 / (kb * param1[0]) - 1 / (kb * param2[0]),
            factor
            * (1 / (kb * param1[0]) * param1[1] - 1 / (kb * param2[0]) * param2[1]),
        ]
    )

    if verbosity > 1:
        print(
            "Analytical slope of {:s}: {:.8f}, {:.8f}".format(
                pstring, trueslope[0], trueslope[1]
            )
        )

    quant = {}

    # ==================== #
    # prepare trajectories #
    # ==================== #
    # Discard burn-in period and time-correlated frames
    traj1 = trajectory.prepare(
        traj1,
        cut=cutoff,
        verbosity=verbosity,
        name="Trajectory 1",
        skip_preparation=data_is_uncorrelated,
    )
    traj2 = trajectory.prepare(
        traj2,
        cut=cutoff,
        verbosity=verbosity,
        name="Trajectory 2",
        skip_preparation=data_is_uncorrelated,
    )

    # calculate overlap
    traj1_full = traj1
    traj2_full = traj2
    traj1, traj2, min_ene, max_ene = trajectory.overlap(
        traj1=traj1_full,
        traj2=traj2_full,
    )
    if verbosity > 0:
        print(
            "Overlap is {:.1%} of trajectory 1 and {:.1%} of trajectory 2.".format(
                traj1.shape[1] / traj1_full.shape[1],
                traj2.shape[1] / traj2_full.shape[1],
            )
        )
    if verbosity > 0:
        cov1 = np.cov(traj1_full)
        sig1 = np.sqrt(np.diag(cov1))
        cov2 = np.cov(traj2_full)
        sig2 = np.sqrt(np.diag(cov2))

        # First parameter is always T
        dt1 = 2 * kb * param1[0] * param1[0] / sig1[0]
        dt2 = 2 * kb * param2[0] * param2[0] / sig2[0]

        # Second parameter is either P or mu
        if dtempdpress:
            parameter_name = "dP"
            # If the second parameter is P, we need to adjust the units
            sig1[1] *= pvconvert
            sig2[1] *= pvconvert
        else:
            parameter_name = "dmu"

        dparam1 = 2 * kb * param1[0] / sig1[1]
        dparam2 = 2 * kb * param2[0] / sig2[1]

        if verbosity > 1:
            print(
                "A rule of thumb states that a good overlap can be expected when choosing state\n"
                "points separated by about 2 standard deviations.\n"
                "For the current trajectories, dT = {:.1f}, and {:s} = {:.1f},\n"
                "with standard deviations sig1 = [{:.1f}, {:.1g}], and sig2 = [{:.1f}, {:.1g}].\n"
                "According to the rule of thumb, given point 1, the estimate is dT = {:.1f}, {:s} = {:.1f}, and\n"
                "                                given point 2, the estimate is dT = {:.1f}, {:s} = {:.1f}.".format(
                    param2[0] - param1[0],
                    parameter_name,
                    abs(param2[1] - param1[1]),
                    sig1[0],
                    sig1[1],
                    sig2[0],
                    sig2[1],
                    dt1,
                    parameter_name,
                    dparam1,
                    dt2,
                    parameter_name,
                    dparam2,
                )
            )

        print(
            "Rule of thumb estimates that (dT,{:s}) = ({:.1f},{:.1f}) would be optimal "
            "(currently, (dT,{:s}) = ({:.1f},{:.1f}))".format(
                parameter_name,
                0.5 * (dt1 + dt2),
                0.5 * (dparam1 + dparam2),
                parameter_name,
                param2[0] - param1[0],
                abs(param2[1] - param1[1]),
            )
        )
    if min_ene is None:
        raise pv_error.InputError(
            ["traj1", "traj2"], "No overlap between trajectories."
        )

    # calculate inefficiency
    g1 = np.array(
        [
            pymbar.timeseries.statisticalInefficiency(traj1[0]),
            pymbar.timeseries.statisticalInefficiency(traj1[1]),
        ]
    )
    g2 = np.array(
        [
            pymbar.timeseries.statisticalInefficiency(traj2[0]),
            pymbar.timeseries.statisticalInefficiency(traj2[1]),
        ]
    )

    w_f = -trueslope[0] * traj1[0] - trueslope[1] * traj1[1]
    w_r = trueslope[0] * traj2[0] + trueslope[1] * traj2[1]

    if verbosity > 2:
        print("Computing log of partition functions using pymbar.BAR...")
    df, ddf = pymbar.BAR(w_f, w_r)
    if verbosity > 2:
        print(
            "Using {:.5f} for log of partition functions as computed from BAR.".format(
                df
            )
        )
        print("Uncertainty in quantity is {:.5f}.".format(ddf))
        print(
            "Assuming this is negligible compared to sampling error at individual points."
        )

    # ================== #
    # max-likelihood fit #
    # ================== #
    if verbosity > 2:
        print("Computing the maximum likelihood parameters")

    fitvals, dfitvals = do_max_likelihood_fit(
        traj1,
        traj2,
        g1,
        g2,
        init_params=np.array([df, trueslope[0], trueslope[1]]),
        verbose=(verbosity > 1),
    )

    slope = fitvals[1:]
    dslope = dfitvals[1:]
    quant["maxLikelihood"] = np.abs((slope - trueslope) / dslope)
    if verbosity > 0:
        print_stats(
            title="Maximum Likelihood Analysis (analytical error)",
            fitvals=fitvals,
            dfitvals=dfitvals,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=None,
            pvconvert=pvconvert,
            dtemp=False,
            dpress=False,
            dmu=False,
            dtempdpress=dtempdpress,
            dtempdmu=dtempdmu,
        )

    if not bootstrap_error:
        return quant["maxLikelihood"]

    # =============================== #
    # bootstrapped max-likelihood fit #
    # =============================== #
    if verbosity > 2:
        print("Computing bootstrapped maximum likelihood parameters")

    if bootstrap_seed is not None:
        np.random.seed(bootstrap_seed)

    bs_fitvals = []
    for t1, t2 in zip(
        trajectory.bootstrap(traj1, bootstrap_repetitions),
        trajectory.bootstrap(traj2, bootstrap_repetitions),
    ):
        # use overlap region
        t1, t2, min_ene, max_ene = trajectory.overlap(traj1=t1, traj2=t2)
        # calculate inefficiency
        g1 = np.array(
            [
                pymbar.timeseries.statisticalInefficiency(t1[0]),
                pymbar.timeseries.statisticalInefficiency(t1[1]),
            ]
        )
        g2 = np.array(
            [
                pymbar.timeseries.statisticalInefficiency(t2[0]),
                pymbar.timeseries.statisticalInefficiency(t2[1]),
            ]
        )
        # calculate max_likelihood fit
        fv, _ = do_max_likelihood_fit(
            t1,
            t2,
            g1,
            g2,
            init_params=np.array([df, trueslope[0], trueslope[1]]),
            verbose=(verbosity > 2),
        )
        bs_fitvals.append(fv)

    bs_fitvals = np.array(bs_fitvals)
    # slope = np.average(fitvals[:, 1:])
    dslope = np.std(bs_fitvals[:, 1:], axis=0)
    quant["bootstrap"] = np.abs((slope - trueslope) / dslope)
    if verbosity > 0:
        print_stats(
            title="Maximum Likelihood Analysis (bootstrapped error)",
            fitvals=np.concatenate(([fitvals], bs_fitvals)),
            dfitvals=None,
            kb=kb,
            param1=param1,
            param2=param2,
            trueslope=trueslope,
            temp=None,
            pvconvert=pvconvert,
            dtemp=False,
            dpress=False,
            dmu=False,
            dtempdpress=dtempdpress,
            dtempdmu=dtempdmu,
        )

    return quant["bootstrap"]
