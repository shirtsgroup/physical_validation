###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (c) 2017-2021 University of Colorado Boulder               #
#              (c) 2012      The University of Virginia                   #
#                                                                         #
###########################################################################

import warnings

import numpy as np
from pymbar import timeseries
from scipy import stats

from . import error as pv_error


def equilibrate(traj, verbose=False, name=None):
    traj = np.array(traj)
    if traj.ndim == 1:
        t0, g, n_eff = timeseries.detectEquilibration(traj)
        if t0 == 0 and traj.size > 10:
            # See https://github.com/choderalab/pymbar/issues/277
            t0x, gx, n_effx = timeseries.detectEquilibration(traj[10:])
            if t0x != 0:
                t0 = t0x + 10
        n = traj.size
        res = traj[t0:]

    elif traj.ndim == 2 and traj.shape[0] == 2:
        t01, g1, n_eff1 = timeseries.detectEquilibration(traj[0])
        t02, g2, n_eff2 = timeseries.detectEquilibration(traj[1])
        t0 = max(t01, t02)
        if t0 == 0 and traj.shape[1] > 10:
            # See https://github.com/choderalab/pymbar/issues/277
            t01x, g1x, n_eff1x = timeseries.detectEquilibration(traj[0, 10:])
            t02x, g2x, n_eff2x = timeseries.detectEquilibration(traj[1, 10:])
            t0x = max(t01x, t02x)
            if t0x != 0:
                t0 = t0x + 10
        n = traj.shape[1]
        res = traj[:, t0:]
    elif traj.ndim == 2:
        raise NotImplementedError(
            "trajectory.equilibrate() in 2 dimensions is only "
            "implemented for exactly two timeseries."
        )
    else:
        raise NotImplementedError(
            "trajectory.equilibrate() is not implemented for "
            "trajectories with more than 2 dimensions."
        )

    if verbose:
        if not name:
            name = "Trajectory"
        if t0 == 0:
            print("{:s} equilibration: No frames discarded for burn-in.".format(name))
        elif t0 == 1:
            print(
                "{:s} equilibration: First frame ({:.1%} of "
                "trajectory) discarded for burn-in.".format(name, 1 / n)
            )
        else:
            print(
                "{:s} equilibration: First {:d} frames ({:.1%} of "
                "trajectory) discarded for burn-in.".format(name, t0, t0 / n)
            )

    return res


def decorrelate(traj, verbose=False, name=None):
    traj = np.array(traj)
    if traj.ndim == 1:
        idx = timeseries.subsampleCorrelatedData(traj)
        n0 = traj.size
        n1 = len(idx)
        res = traj[idx]
    elif traj.ndim == 2:
        # pymbar doesn't offer to decorrelate two samples, so let's do it ourselves
        # and just use the decorrelation of the sample more strongly correlated
        #
        # calculate (maximal) inefficiency
        g1 = timeseries.statisticalInefficiency(traj[0])
        g2 = timeseries.statisticalInefficiency(traj[1])
        g = np.max([g1, g2])
        # calculate index
        n0 = traj.shape[1]
        idx = np.unique(
            np.array(np.round(np.arange(0, int(n0 / g + 0.5)) * g), dtype=int)
        )
        idx = idx[idx < n0]
        n1 = len(idx)
        res = traj[:, idx]
    else:
        raise NotImplementedError(
            "trajectory.decorrelate() is not implemented for "
            "trajectories with more than 1 dimension."
        )
    if verbose:
        n = n0 - n1
        if not name:
            name = "Trajectory"
        if n == 0:
            print(
                "{:s} decorrelation: No frames discarded for decorrelation.".format(
                    name
                )
            )
        elif n == 1:
            print(
                "{:s} decorrelation: 1 frame ({:.1%} of "
                "trajectory) discarded for decorrelation.".format(name, 1 / n0)
            )
        else:
            print(
                "{:s} decorrelation: {:d} frames ({:.1%} of "
                "trajectory) discarded for decorrelation.".format(name, n, n / n0)
            )

    return res


def cut_tails(traj, cut, verbose=False, name=None):
    traj = np.array(traj)
    dc = 100 * cut
    if traj.ndim == 1:
        with warnings.catch_warnings():
            # With some combination of python version / scipy version,
            # scoreatpercentile throws a warning
            warnings.filterwarnings("ignore", category=FutureWarning)
            tmax = stats.scoreatpercentile(traj, 100 - dc)
            tmin = stats.scoreatpercentile(traj, dc)
        t = traj[(tmin <= traj) * (traj <= tmax)]
        n0 = traj.size
        n = t.size
    elif traj.ndim == 2:
        with warnings.catch_warnings():
            # With some combination of python version / scipy version,
            # scoreatpercentile throws a warning
            warnings.filterwarnings("ignore", category=FutureWarning)
            tmax = stats.scoreatpercentile(traj, 100 - dc, axis=1)
            tmin = stats.scoreatpercentile(traj, dc, axis=1)
        t = traj[
            :,
            (tmin[0] <= traj[0])
            * (tmin[1] <= traj[1])
            * (tmax[0] >= traj[0])
            * (tmax[1] >= traj[1]),
        ]
        n0 = traj.shape[1]
        n = t.shape[1]
    else:
        raise NotImplementedError(
            "trajectory.cut_tails() is not implemented for "
            "trajectories with more than 2 dimension."
        )

    if verbose:
        if not name:
            name = "Trajectory"
        print(
            "{:s} tails (cut = {:.2%}): {:n} frames ({:.2%} of trajectory) were cut".format(
                name, cut, n0 - n, (n0 - n) / n0
            )
        )

    return t


def prepare(traj, cut=None, verbosity=1, name=None):
    traj = np.array(traj)
    if not name:
        name = "Trajectory"

    def traj_length(t):
        if t.ndim == 1:
            return t.size
        else:
            return t.shape[1]

    if traj.ndim > 2:
        raise NotImplementedError(
            "trajectory.prepare() is not implemented for "
            "trajectories with more than 2 dimensions."
        )

    # original length
    n0 = traj_length(traj)
    # equilibrate
    res = equilibrate(traj, verbose=False)
    n1 = traj_length(res)
    if verbosity > 2:
        print(
            "{:s} equilibration: First {:d} frames ({:.1%} of "
            "trajectory) discarded for burn-in.".format(name, n0 - n1, (n0 - n1) / n0)
        )
    # decorrelate
    res = decorrelate(res, verbose=False)
    n2 = traj_length(res)
    if verbosity > 2:
        print(
            "{:s} decorrelation: {:d} frames ({:.1%} of equilibrated "
            "trajectory) discarded for decorrelation.".format(
                name, n1 - n2, (n1 - n2) / n1
            )
        )
    # cut tails
    if cut is not None:
        res = cut_tails(res, cut, verbose=False)
        n3 = traj_length(res)
        if verbosity > 2:
            print(
                "{:s} tails (cut = {:.2%}): {:n} frames ({:.2%} of equilibrated and "
                "decorrelated trajectory) were cut".format(
                    name, cut, n2 - n3, (n2 - n3) / n2
                )
            )
    # end length
    nn = traj_length(res)

    if verbosity > 0:
        print(
            "After equilibration, decorrelation and tail pruning, {:.2%} ({:n} frames) "
            "of original {:s} remain.".format(nn / n0, nn, name)
        )

    return res


def overlap(traj1, traj2, cut=None, verbose=False, name=None):
    traj1 = np.array(traj1)
    traj2 = np.array(traj2)
    if traj1.ndim == traj2.ndim and traj2.ndim == 1:
        if cut:
            dc = 100 * cut
            max1 = stats.scoreatpercentile(traj1, 100 - dc)
            min1 = stats.scoreatpercentile(traj1, dc)
            max2 = stats.scoreatpercentile(traj2, 100 - dc)
            min2 = stats.scoreatpercentile(traj2, dc)
        else:
            max1 = traj1.max()
            min1 = traj1.min()
            max2 = traj2.max()
            min2 = traj2.min()

        tmin = max(min1, min2)
        tmax = min(max1, max2)

        t1 = traj1[(tmin <= traj1) * (traj1 <= tmax)]
        t2 = traj2[(tmin <= traj2) * (traj2 <= tmax)]
    elif traj1.ndim == traj2.ndim and traj2.ndim == 2:
        if traj1.shape[0] != 2 or traj2.shape[0] != 2:
            raise NotImplementedError(
                "trajectory.overlap() in 2 dimensions is only "
                "implemented for exactly two timeseries per trajectory."
            )
        if cut:
            dc = 100 * cut
            max1 = stats.scoreatpercentile(traj1, 100 - dc, axis=1)
            min1 = stats.scoreatpercentile(traj1, dc, axis=1)
            max2 = stats.scoreatpercentile(traj2, 100 - dc, axis=1)
            min2 = stats.scoreatpercentile(traj2, dc, axis=1)
        else:
            max1 = traj1.max(axis=1)
            min1 = traj1.min(axis=1)
            max2 = traj2.max(axis=1)
            min2 = traj2.min(axis=1)

        tmin = np.max([min1, min2], axis=0)
        tmax = np.min([max1, max2], axis=0)

        t1 = traj1[
            :,
            (tmin[0] <= traj1[0])
            * (tmin[1] <= traj1[1])
            * (tmax[0] >= traj1[0])
            * (tmax[1] >= traj1[1]),
        ]
        t2 = traj2[
            :,
            (tmin[0] <= traj2[0])
            * (tmin[1] <= traj2[1])
            * (tmax[0] >= traj2[0])
            * (tmax[1] >= traj2[1]),
        ]
    elif traj1.ndim != traj2.ndim:
        raise pv_error.InputError(
            ["traj1", "traj2"], "Trajectories don't have the same number of dimensions"
        )
    else:
        raise NotImplementedError(
            "trajectory.overlap() is not implemented for "
            "trajectories with more than 2 dimensions."
        )

    if np.any(max1 < min2) or np.any(max2 < min1):
        if verbose:
            if not name:
                name = "Trajectory"
            print("{:s} overlap: No overlap found between trajectories".format(name))
        return np.array([]), np.array([]), None, None

    if verbose:
        if not name:
            name = "Trajectory"
        print(
            "{:s} overlap: {:.1%} of trajectory 1, and {:.1%} of trajectory 2 "
            "were found within overlap region.\n"
            "              That corresponds to {:n} frames and {:n} frames, "
            "respectively".format(
                name, len(traj1) / len(t1), len(traj2) / len(t2), len(t1), len(t2)
            )
        )

    return t1, t2, tmin, tmax


def bootstrap(traj, n_samples):
    traj = np.array(traj)
    if traj.ndim == 1:
        n_traj = traj.size
    elif traj.ndim == 2:
        n_traj = traj.shape[1]
    else:
        raise NotImplementedError(
            "trajectory.bootstrap() is not implemented for "
            "trajectories with more than 2 dimensions."
        )
    for _ in range(n_samples):
        resample_idx = np.floor(np.random.rand(n_traj) * n_traj).astype(int)
        if traj.ndim == 1:
            yield traj[resample_idx]
        else:
            yield traj[:, resample_idx]


def jackknife(traj):
    traj = np.array(traj)
    if traj.ndim == 1:
        n_traj = traj.size
        for idx in range(n_traj):
            yield np.concatenate((traj[:idx], traj[idx + 1 :]))
    elif traj.ndim == 2:
        n_traj = traj.shape[1]
        for idx in range(n_traj):
            yield np.concatenate((traj[:, :idx], traj[:, idx + 1 :]))
    else:
        raise NotImplementedError(
            "trajectory.jackknife() is not implemented for "
            "trajectories with more than 2 dimensions."
        )


def bca(param, param_boot, param_jack, alpha):
    param_boot = np.array(param_boot)
    param_jack = np.array(param_jack)
    n_boot = param_boot.size

    z0 = stats.norm.ppf(np.sum(param_boot < param) / n_boot)
    zlow = stats.norm.ppf(alpha / 2)
    zhigh = stats.norm.ppf(1 - alpha / 2)

    avg_p_jack = np.mean(param_jack)
    a = np.sum((avg_p_jack - param_jack) ** 3) / (
        6 * (np.sum((avg_p_jack - param_jack) ** 2) ** (3 / 2))
    )

    a1 = stats.norm.cdf(z0 + (z0 + zlow) / (1 - a * (z0 + zlow)))
    a2 = stats.norm.cdf(z0 + (z0 + zhigh) / (1 - a * (z0 + zhigh)))

    return np.percentile(param_boot, a1 * 100), np.percentile(param_boot, a2 * 100)
