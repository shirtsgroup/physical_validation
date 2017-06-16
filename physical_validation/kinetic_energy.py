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
"""
The `kinetic_energy` module is part of the physical_validation package, and
consists of checks of the kinetic energy distribution and its
equipartition.
"""
from __future__ import print_function
from __future__ import division

import physical_validation.util.kinetic_energy as util_kin


def mb_ensemble(data, alpha, verbose=False):
    r"""Checks if a kinetic energy trajectory is Maxwell-Boltzmann distributed.

    Parameters
    ----------
    data : SimulationData
        Simulation data object
    alpha : float
        Confidence. TODO: Check proper statistical definition.
    verbose : bool, optional
        Print result details. Default: False.
    
    Returns
    -------
    result : bool
        True if hypothesis stands, false if hypothesis got rejected.
        
    Notes
    -----
    This function checks whether the hypothesis that a sample 
    of kinetic energies is Maxwell-Boltzmann distributed given a specific
    target temperature and the number of degrees of freedom in the system,

    .. math::
        P(K) \sim K^{N/2} e^{-\beta K} \, ,

    holds under a given confidence level :math:`\alpha`.
    The check is performed using the Kolmogorov-Smirnov test provided by
    scipy.stats.kstest_.

    .. _scipy.stats.kstest: https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.kstest.html
    
    .. note:: The Kolmogorov-Smirnov is known to have two weaknesses.
    
       #. The test is more sensitive towards deviations around the center
          of the distribution than at its tails. We deem this to be acceptable
          for most MD applications, but be wary if yours is sensible to the
          kinetic distribution tails.
       #. The test is not valid if its parameters are guessed from the data
          set. Using the target temperature of the MD simulation as an input 
          is therefore perfectly valid, but using the average temperature
          over the trajectory as an input to the test can potentially 
          invalidate it.

    .. todo:: Can we check the influence of sample size on test results?

    """
    return util_kin.check_mb_ensemble(data.observables['kinetic_energy'],
                                      data.ensemble.temperature,
                                      data.topology.ndof_total, alpha,
                                      data.units.kb, verbose)


def equipartition(data, dtemp=0.1, temp=None, alpha=0.05,
                  molec_groups=None,
                  random_divisions=0, random_groups=0,
                  verbosity=2):
    r"""Checks the equipartition of a simulation trajectory.
    
    Parameters
    ----------
    data : SimulationData
        Simulation data object
    dtemp : float, optional
        Fraction of temperature deviation tolerated between groups. Default: 0.05 (5%).
    temp : float, optional
        Target temperature of the simulation. If None, the kinetic energies will not be
        tested for Maxwell-Boltzmann distribution, but only compared amongst each others.
        Default: None.
    alpha : float, optional
        Confidence for Maxwell-Boltzmann test. Default: 0.05 (5%).
    molec_groups : list of array-like (ngroups x ?), optional
        List of 1d arrays containing molecule indeces defining groups. Useful to pre-define
        groups of molecules (e.g. solute / solvent, liquid mixture species, ...). If None,
        no pre-defined molecule groups will be tested. Default: None.
    random_divisions : int, optional
        Number of random division tests attempted. Default: 0 (random division tests off).
    random_groups : int, optional
        Number of groups the system is randomly divided in. Default: 2.
    verbosity : int, optional
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 2.

    Returns
    -------
    result : int
        Number of equipartition violations. Tune up verbosity for details.
        
    Notes
    -----
    
    
    This function compares the kinetic energy between groups of degrees of
    freedom. Theoretically, the kinetic energy is expected (via the 
    equipartition theorem) to be equally distributed over all degrees of 
    freedom. In practice, deviations of temperature between groups of 
    degrees of freedom up to several degrees K are routinely observed. 
    Larger deviations can, however, hint to misbehaving simulations, such 
    as, e.g., frozen degrees of freedom, lack of energy exchange between 
    degrees of freedom, and transfer of heat from faster to slower 
    oscillating degrees of freedom.
    
    Splitting of degrees of freedom is done both on a sub-molecular and on
    a molecular level. On a sub-molecular level, the degrees of freedom of
    a molecule can be partitioned into rigid-body contributions 
    (translation of the center-of-mass, rotation around the 
    center-of-mass) and intra-molecular contributions. On a molecular 
    level, the single molecules of the system can be divided in groups, 
    either by function (solute / solvent, different species of liquid 
    mixtures, ...) or randomly.
    
    `check_equipartition()` compares the partitioned temperatures of the 
    entire system and, optionally, of predefined or randomly separated 
    groups. 
    
    Note: In theory, the kinetic energy of the subgroups are expected to 
    be individually Maxwell-Boltzmann distributed. As this is seldomly 
    holding in practice (see above), `check_equipartition()` is by 
    default checking only for abnormal deviations in average temperatures.
    The more strict Maxwell-Boltzmann testing can be invoked by giving the
    target temperature `temp` as an input.

    """
    return util_kin.check_equipartition(positions=data.trajectory['position'],
                                        velocities=data.trajectory['velocity'],
                                        masses=data.topology.mass,
                                        molec_idx=data.topology.molecule_idx,
                                        molec_nbonds=data.topology.nconstraints_per_molecule,
                                        ndof_reduction_tra=data.topology.ndof_reduction_tra,
                                        ndof_reduction_rot=data.topology.ndof_reduction_rot,
                                        dtemp=dtemp, temp=temp, alpha=alpha,
                                        molec_groups=molec_groups,
                                        random_divisions=random_divisions,
                                        random_groups=random_groups,
                                        verbosity=verbosity)
