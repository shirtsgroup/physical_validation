"""
integratorconvergence.py
@author Pascal T. Merz
"""
# py2 compatibility
from __future__ import absolute_import, division, print_function
# other imports
import numpy as np

import physicalvalidation.util.integrator as util_integ
from physicalvalidation.data.simulation_data import SimulationData
import physicalvalidation.util.error as pv_error


def convergence(simulations,
                convergence_test=util_integ.simple_convergence_test,
                verbose=True, slope=False, tol=0.1):
    r"""
    
    Parameters
    ----------
    simulations : list of SimulationData
    convergence_test : callable, optional
    verbose : bool, optional
    slope : bool, optional
    tol : float, optional

    Returns
    -------
    result : bool

    """
    constant_of_motion = {}

    for s in simulations:
        if not isinstance(s, SimulationData):
            raise pv_error.InputError('simulations',
                                      'Expected a list of objects of type SimulationData')
        if s.dt <= 0:
            raise pv_error.InputError('simulations',
                                      'Found SimulationData with timestep dt=0')
        key = str(s.dt)

        if key in constant_of_motion:
            raise pv_error.InputError('simulations',
                                      'Found two simulations with identical timestep')

        constant_of_motion[key] = s.observables.constant_of_motion

    return util_integ.check_convergence(constant_of_motion,
                                        convergence_test=convergence_test,
                                        verbose=verbose,
                                        slope=slope,
                                        tol=tol)
