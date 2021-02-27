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
r"""
This module contains low-level functionality of the
`physical_validation.integrator` module. The functions in this module should
generally not be called directly. Please use the high-level functions from
`physical_validation.integrator`.
"""
import numpy as np

from . import plot


def calculate_rmsd(data, time=None, slope=False):
    assert isinstance(data, np.ndarray) and data.ndim == 1
    assert time is None or isinstance(time, np.ndarray) and time.ndim == 1

    avg = data.mean()

    if time is None:
        time = np.arange(data.size)

    fit = np.polyfit(time, data, 1)

    def f(x):
        return fit[0] * x + fit[1]

    if slope:
        rmsd = 0
        for t, d in zip(time, data):
            rmsd += (d - f(t)) ** 2
        rmsd = np.sqrt(rmsd / data.size)
    else:
        rmsd = data.std()

    return avg, rmsd, fit[0]


def max_deviation(dts, rmsds):
    dt_ratio_2 = (dts[:-1] / dts[1:]) ** 2
    rmsds = rmsds[:-1] / rmsds[1:]
    return np.max(np.abs(1 - rmsds / dt_ratio_2))


def check_convergence(
    const_traj,
    convergence_test=max_deviation,
    verbose=True,
    slope=False,
    screen=False,
    filename=None,
):

    assert isinstance(const_traj, dict)
    assert len(const_traj) >= 2

    if verbose:
        print("{:65s}".format("-" * 65))
        print(
            "{:>10s} {:>10s} {:>10s} {:>10s} {:^21s}".format(
                "dt", "avg", "rmsd", "slope", "ratio"
            )
        )
        print("{:43s} {:>10s} {:>10s}".format("", "dt^2", "rmsd"))
        print("{:65s}".format("-" * 65))

    prev = None

    results = {}
    for dt, traj in sorted(const_traj.items(), key=lambda x: float(x[0]), reverse=True):
        assert isinstance(traj, np.ndarray)
        assert traj.ndim == 1 or traj.ndim == 2
        dt = float(dt)

        if traj.ndim == 1:
            data = traj
            time = None
        else:
            data = traj[1]
            time = traj[0]

        results[dt] = calculate_rmsd(data, time, slope)

        if verbose:
            if prev is None:
                print(
                    "{:10.4g} {:10.2f} {:10.2e} {:10.2e} {:>10s} {:>10s}".format(
                        dt, results[dt][0], results[dt][1], results[dt][2], "--", "--"
                    )
                )
                prev = [dt, results[dt][1]]
            else:
                print(
                    "{:10.4g} {:10.2f} {:10.2e} {:10.2e} {:10.2f} {:10.2f}".format(
                        dt,
                        results[dt][0],
                        results[dt][1],
                        results[dt][2],
                        prev[0] ** 2 / dt ** 2,
                        prev[1] / results[dt][1],
                    )
                )
                prev = [dt, results[dt][1]]

    if verbose:
        print("{:65s}".format("-" * 65))

    dts = np.sort(np.array([float(dt) for dt in results.keys()]))[::-1]
    rmsds = np.array([float(results[dt][1]) for dt in dts])

    do_plot = screen or filename is not None

    if do_plot:
        data = [
            {
                "x": dts[1:],
                "y": rmsds[:-1] / rmsds[1:],
                "name": "Integrator convergence",
            },
            {
                "x": dts[1:],
                "y": (dts[:-1] / dts[1:]) ** 2,
                "name": "Expected convergence",
            },
        ]

        plot.plot(
            data,
            legend="best",
            title="Actual vs. expected convergence",
            xlabel="Time step",
            ylabel="Convergence",
            xlim=(0, dts[1]),
            inv_x=True,
            filename=filename,
            screen=screen,
        )

    return convergence_test(dts, rmsds)
