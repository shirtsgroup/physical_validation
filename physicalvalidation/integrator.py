###########################################################################
#                                                                         #
#    physicalvalidation,                                                  #
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
"""
The `integratorconvergence` module is part of the physicalvalidation
package, and consists of checks of the convergence of the MD integrator.
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
