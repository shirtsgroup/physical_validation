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
"""
The `integrator.convergence` module is part of the physical_validation
package, and consists of checks of the convergence of the MD integrator.
"""
from typing import List, Optional

from .data import SimulationData
from .util import error as pv_error
from .util import integrator as util_integ


def convergence(
    simulations: List[SimulationData],
    convergence_test: str = "max_deviation",
    verbose: bool = True,
    screen: bool = False,
    filename: Optional[str] = None,
) -> float:
    r"""
    Compares the convergence of the fluctuations of conserved quantities
    with decreasing simulation time step to theoretical expectations.

    Parameters
    ----------
    simulations
        The (otherwise identical) simulations performed using different
        time steps
    convergence_test
        A function defining the convergence test. Currently, only one
        test is implemented:
        `max_deviation`, which is chosen by default
    verbose
        If True, print more detailed output. Default: False.
    screen
        Plot convergence on screen. Default: False.
    filename
        Plot convergence to `filename`. Default: None, no plotting to file.

    Returns
    -------
        The largest deviation from the expected ratio of fluctuations.

    Notes
    -----
    For a simplectic integration algorithm, the fluctuations
    :math:`\delta E` of a constant of motion :math:`E` (such as, for
    example, the total energy in a NVE simulations) are theoretically
    expected to scale like the squared timestep of the integration.
    When comparing two otherwise identical simulations performed at
    different time step :math:`\Delta t`, the following equality is
    hence expected to hold:

    .. math::
        \frac{\Delta t_1^{2}}{\Delta t_2^{2}} = \frac{\delta E_1}{\delta E_2}

    This function calculates the ratio of the fluctuation for simulations
    performed at different timesteps and compares it to the analytically
    expected value. If the deviation is larger than `tol`, the test is
    considered failed.

    """
    constant_of_motion = {}

    convergence_tests = {"max_deviation": util_integ.max_deviation}

    if convergence_test not in convergence_tests:
        raise pv_error.InputError("convergence_test", "Unknown convergence test.")

    convergence_test = convergence_tests[convergence_test]

    for s in simulations:
        if not isinstance(s, SimulationData):
            raise pv_error.InputError(
                "simulations", "Expected a list of objects of type SimulationData"
            )
        if s.dt <= 0:
            raise pv_error.InputError(
                "simulations", "Found SimulationData with timestep dt<=0"
            )
        key = str(s.dt)

        if key in constant_of_motion:
            raise pv_error.InputError(
                "simulations", "Found two simulations with identical timestep"
            )

        s.raise_if_observable_data_is_invalid(
            required_observables=["constant_of_motion"],
            test_name="integrator.convergence",
            argument_name="simulations",
        )
        constant_of_motion[key] = s.observables.constant_of_motion

    return util_integ.check_convergence(
        constant_of_motion,
        convergence_test=convergence_test,
        verbose=verbose,
        screen=screen,
        filename=filename,
    )
