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
r"""
This module contains low-level functionality of the
`physicalvalidation.kineticenergy` module. The functions in this module
should generally not be called directly. Please use the high-level
functions from `physicalvalidation.kinetic energy`.
"""
from __future__ import print_function
from __future__ import division

import scipy.stats as stats
import numpy as np


def temperature(kin, ndof, kb=8.314e-3):
    r"""
    Calculates the temperature acccording to the equipartition theorem.
    
    .. math::
        T(K) = \frac{2K}{N k_B}

    Parameters
    ----------
    kin : float
        Kinetic energy.
    ndof : float
        Number of degrees of freedom.
    kb : float
        Boltzmann constant :math:`k_B`. Default: 8.314e-3 (kJ/mol).
    
    Returns
    -------
    temperature : float
        Calculated temperature.
    """
    # ndof * kb * T = 2 * kin
    return 2 * float(kin) / (float(ndof) * float(kb))


def check_mb_ensemble(kin, temp, ndof, alpha, kb=8.314e-3, verbose=False,
                      plot=False, plot_name=None):
    r"""
    Checks if a kinetic energy trajectory is Maxwell-Boltzmann distributed.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using 
       the high-level version based on the SimulationData object. See
       physicalvalidation.kineticenergy.check_mb_ensemble for more
       information and full documentation.
       
    Parameters
    ----------
    kin : array-like
        Kinetic energy snapshots of the system.
    temp : float
        Target temperature of the system. Used to construct the 
        Maxwell-Boltzmann distribution.         
    ndof : float
        Number of degrees of freedom in the system. Used to construct the 
        Maxwell-Boltzmann distribution.
    alpha : float
        Confidence. TODO: Check proper statistical definition.
    kb : float
        Boltzmann constant :math:`k_B`. Default: 8.314e-3 (kJ/mol).
    verbose : bool
        Print result details. Default: False.
    plot : bool
        Only for debug - prints a plottable file.
    plot_name : string
        Only for debug - the name of the file to write the plotting data to.
        
    Returns
    -------
    result : bool
        True if hypothesis stands, false if hypothesis got rejected.
        
    See Also
    --------
    physicalvalidation.kineticenergy.check_mb_ensemble : High-level version
    """

    kt = kb * temp
    d, p = stats.kstest(kin, 'chi2', (ndof, 0, kt/2))

    if plot:
        h, k = np.histogram(kin, bins=25, normed=True)
        k = (k[1:] + k[:-1])/2
        with open(plot_name + '.dat', 'w') as f:
            for hh, kk in zip(h, k):
                f.write('{:f} {:f}\n'.format(kk, hh))
        ana = stats.chi2(df=ndof, scale=kt/2)
        k = np.linspace(ana.ppf(0.0001), ana.ppf(0.9999), 100)
        with open(plot_name + '_ana.dat', 'w') as f:
            for kk in k:
                f.write('{:f} {:f}\n'.format(kk, ana.pdf(kk)))

    if p > alpha:
        if verbose:
            print('Kolmogorov-Smirnov test result: p = {:f}\n'
                  'Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed\n'
                  'Confidence alpha = {:f}\n'
                  'Result: Hypothesis stands'.format(p, alpha))
        return True

    if verbose:
        print('Kolmogorov-Smirnov test result: p = {:f}\n'
              'Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed\n'
              'Confidence alpha = {:f}\n'
              'Result: Hypothesis rejected'.format(p, alpha))
    return False


def check_equipartition(positions, velocities, masses,
                        molec_idx, molec_nbonds,
                        natoms, nmolecs, nframes,
                        ndof_reduction_tra=0, ndof_reduction_rot=0,
                        dtemp=0.1, temp=None, alpha=0.05,
                        molec_groups=None,
                        random_divisions=0, random_groups=2,
                        verbosity=2):
    r"""
    Checks the equipartition of a simulation trajectory.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using 
       the high-level version based on the SimulationData object. See
       physicalvalidation.kineticenergy.check_equipartition for more
       information and full documentation.

    Parameters
    ----------    
    positions : array-like (nframes x natoms x 3)
        3d array containing the positions of all atoms for all frames
    velocities : array-like (nframes x natoms x 3)
        3d array containing the velocities of all atoms for all frames
    masses : array-like (natoms x 1)
        1d array containing the masses of all atoms
    molec_idx : array-like (nmolecs x 1)
        Index of first atom for every molecule 
    molec_nbonds : array-like (nmolecs x 1)
        Number of bonds for every molecule
    natoms : int
        Number of atoms in the system
    nmolecs : int
        Number of molecules in the system
    nframes : int
        Number of frames contained in the `positions` and `velocities`
        vectors
    ndof_reduction_tra : int, optional
        Number of center-of-mass translational degrees of freedom to 
        remove. Default: 0.
    ndof_reduction_rot : int, optional
        Number of center-of-mass rotational degrees of freedom to remove. 
        Default: 0.
    dtemp : float, optional
        Fraction of temperature deviation tolerated between groups. 
        Default: 0.05 (5%).
    temp : float, optional
        Target temperature of the simulation. If None, the kinetic 
        energies will not be tested for Maxwell-Boltzmann distribution, 
        but only compared amongst each others. Default: None.
    alpha : float, optional
        Confidence for Maxwell-Boltzmann test. Default: 0.05 (5%).
    molec_groups : list of array-like (ngroups x ?), optional
        List of 1d arrays containing molecule indeces defining groups. 
        Useful to pre-define groups of molecules (e.g. solute / solvent, 
        liquid mixture species, ...). If None, no pre-defined molecule 
        groups will be tested. Default: None.
    random_divisions : int, optional
        Number of random division tests attempted. Default: 0 (random 
        division tests off).
    random_groups : int, optional
        Number of groups the system is randomly divided in. Default: 2.
    verbosity : int, optional
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 2.

    Returns
    -------
    result : int
        Number of equipartition violations. Tune up verbosity for details.
        
    See Also
    --------
    physicalvalidation.kineticenergy.check_equipartition : High-level version

    """

    dict_keys = ['tot', 'tra', 'rni', 'rot', 'int']

    # for each molecule, calculate total / translational / rotational & internal /
    #   rotational / internal degrees of freedom
    #   returns: list of dict of floats (shape: nmolecs x 5 x 1)
    ndof_molec = calc_ndof(natoms, nmolecs, molec_idx, molec_nbonds,
                           ndof_reduction_tra, ndof_reduction_rot)

    # for each frame, calculate total / translational / rotational & internal /
    #   rotational / internal kinetic energy for each molecule
    kin_molec = []
    for r, v in zip(positions, velocities):
        kin_molec.append(calc_kinetic_energy(r, v, masses,
                                             molec_idx, natoms, nmolecs))

    result = 0

    # test system-wide tot, tra, rni, rot, int
    if temp is not None:
        # check for Maxwell-Boltzmann distribution of
        # partitioned kinetic energy trajectories
        result += test_mb_dist(kin_molec, ndof_molec, nmolecs,
                               temp, alpha, dict_keys, verbosity)
    else:
        # compare partitioned temperatures to total temperature
        result += test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                 dtemp, dict_keys, verbosity)

    # divide in random groups
    for i in range(random_divisions):
        # randomly assign a group index to each molecule
        group_idx = np.random.randint(random_groups, size=nmolecs)
        # create molecule index for each group
        groups = []
        for rg in range(random_groups):
            groups.append(np.arange(nmolecs)[group_idx == rg])
        # test each group separately
        for rg, group in enumerate(groups):
            if verbosity > 0:
                print('Testing randomly divided group {:d}'.format(rg))
            if verbosity > 3:
                print(group)
            if temp is not None:
                result += test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                       alpha, dict_keys, group, verbosity)
            else:
                result += test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                         dtemp, dict_keys, group, verbosity)
        # test groups against each others
        for rg1, group1 in enumerate(groups):
            for group2 in groups[rg1 + 1:]:
                test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                      group1, group2,
                                      dtemp, dict_keys, verbosity)

    # use predefined group division
    if molec_groups is not None:
        for mg, group in enumerate(molec_groups):
            if verbosity > 0:
                print('Testing predifined divided group {:d}'.format(mg))
            if verbosity > 3:
                print(group)
            if temp is not None:
                result += test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                       alpha, dict_keys, group, verbosity)
            else:
                result += test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                         dtemp, dict_keys, group, verbosity)
        # test groups against each others
        for rg1, group1 in enumerate(molec_groups):
            for group2 in molec_groups[rg1 + 1:]:
                test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                      group1, group2,
                                      dtemp, dict_keys, verbosity)

    return result


def calc_system_ndof(natoms, nmolecs, nbonds,
                     stop_com_tra, stop_com_rot):
    r"""
    Calculates the total / translational / rotational & internal / 
    rotational / internal degrees of freedom of the system.

    Parameters
    ----------
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system
    nbonds : int
        Total number of bonds in the system
    stop_com_tra : bool
        Was the center-of-mass translation removed during the simulation?
    stop_com_rot : bool
        Was the center-of-mass translation removed during the simulation?
    
    Returns
    -------
    ndof : dict
        Dictionary containing the degrees of freedom.
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # total ndof
    ndof_tot = 3*natoms - nbonds

    # ndof reduction due to COM motion constraining
    if stop_com_tra:
        ndof_tot -= 3
    if stop_com_rot:
        ndof_tot -= 3

    # translational ndof
    ndof_tot_tra = 3*nmolecs
    if stop_com_tra:
        ndof_tot -= 3

    # rotational & internal ndof
    ndof_tot_rni = ndof_tot - ndof_tot_tra

    # rotational ndof
    ndof_tot_rot = 3*nmolecs
    if stop_com_tra:
        ndof_tot -= 3

    # internal ndof
    ndof_tot_int = ndof_tot_rni - ndof_tot_rot

    # return dict
    ndof = {'tot': ndof_tot,
            'tra': ndof_tot_tra,
            'rni': ndof_tot_rni,
            'rot': ndof_tot_rot,
            'int': ndof_tot_int}
    return ndof


def calc_ndof(natoms, nmolecs,
              molec_idx, molec_nbonds,
              ndof_reduction_tra, ndof_reduction_rot):
    r"""
    Calculates the total / translational / rotational & internal / 
    rotational / internal degrees of freedom per molecule.

    Parameters
    ----------
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system
    molec_idx : list of int
        Index of first atom for every molecule 
    molec_nbonds : list of int
        Number of bonds for every molecule
    ndof_reduction_tra : int\
        Number of center-of-mass translational degrees of freedom to 
        remove. Default: 0.
    ndof_reduction_rot : int
        Number of center-of-mass rotational degrees of freedom to remove. 
        Default: 0.

    Returns
    -------
    ndof_molec : list of dict
        List of dictionaries containing the degrees of freedom for each molecule
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # ndof to be deducted per molecule
    # ndof reduction due to COM motion constraining
    ndof_com_tra_pm = ndof_reduction_tra / nmolecs
    ndof_com_rot_pm = ndof_reduction_rot / nmolecs

    ndof_molec = []
    # add last idx to molec_idx to ease looping
    molec_idx = np.append(molec_idx, [natoms])
    # loop over molecules
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(zip(molec_idx[:-1], molec_idx[1:])):
        natoms = idx_atm_end - idx_atm_init
        nbonds = molec_nbonds[idx_molec]
        ndof_tot = 3*natoms - nbonds - ndof_com_tra_pm - ndof_com_tra_pm
        ndof_tra = 3 - ndof_com_tra_pm
        ndof_rot = 3 - ndof_com_rot_pm
        ndof_molec.append({'tot': ndof_tot,
                           'tra': ndof_tra,
                           'rni': ndof_tot - ndof_tra,
                           'rot': ndof_rot,
                           'int': ndof_tot - ndof_tra - ndof_rot})

    return ndof_molec


def calc_kinetic_energy(pos, vel, masses,
                        molec_idx, natoms, nmolecs):
    r"""
    Calculates the total / translational / rotational & internal / 
    rotational / internal kinetic energy per molecule. 

    Parameters
    ----------
    pos : nd-array (natoms x 3)
        2d array containing the positions of all atoms
    vel : nd-array (natoms x 3)
        2d array containing the velocities of all atoms
    masses : nd-array (natoms x 1)
        1d array containing the masses of all atoms
    molec_idx : nd-array (nmolecs x 1)
        Index of first atom for every molecule 
    natoms : int
        Total number of atoms in the system
    nmolecs : int
        Total number of molecules in the system
    
    Returns
    -------
    kin : list of dict
        List of dictionaries containing the kinetic energies for each molecule
        Keys: ['tot', 'tra', 'rni', 'rot', 'int']
    """
    # add last idx to molec_idx to ease looping
    molec_idx = np.append(molec_idx, [natoms])

    # calculate kinetic energy
    kin_tot = np.zeros(nmolecs)
    kin_tra = np.zeros(nmolecs)
    kin_rni = np.zeros(nmolecs)
    kin_rot = np.zeros(nmolecs)
    kin_int = np.zeros(nmolecs)
    # loop over molecules
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(zip(molec_idx[:-1], molec_idx[1:])):
        # compute center of mass position, velocity and total mass
        com_r = np.zeros(3)
        com_v = np.zeros(3)
        com_m = 0
        # loop over atoms in molecule
        for r, v, m in zip(pos[idx_atm_init:idx_atm_end],
                           vel[idx_atm_init:idx_atm_end],
                           masses[idx_atm_init:idx_atm_end]):
            com_r += m*r
            com_v += m*v
            com_m += m

            # total kinetic energy is straightforward
            kin_tot[idx_molec] += .5 * m * np.dot(v, v)

        com_r /= com_m
        com_v /= com_m

        # translational kinetic energy
        kin_tra[idx_molec] = .5 * com_m * np.dot(com_v, com_v)
        # combined rotational and internal kinetic energy
        kin_rni[idx_molec] = kin_tot[idx_molec] - kin_tra[idx_molec]

        # compute tensor of inertia and angular momentum
        inertia = np.zeros((3, 3))
        angular_mom = np.zeros(3)
        # loop over atoms in molecule
        for r, v, m in zip(pos[idx_atm_init:idx_atm_end],
                           vel[idx_atm_init:idx_atm_end],
                           masses[idx_atm_init:idx_atm_end]):
            # relative positions and velocities
            r -= com_r
            v -= com_v
            r2 = np.dot(r, r)
            # inertia tensor:
            #   (i,i) = m*(r*r - r(i)*r(i))
            #   (i,j) = m*r(i)*r(j) (i != j)
            atm_inertia = -m*np.tensordot(r, r, axes=0)
            for i in range(3):
                atm_inertia[i][i] += m*r2
            inertia += atm_inertia
            # angular momentum: r x p
            angular_mom += m * np.cross(r, v)

        # angular velocity of the molecule: inertia^{-1} * angular_mom
        angular_v = np.dot(np.linalg.inv(inertia), angular_mom)

        kin_rot[idx_molec] = .5 * np.dot(angular_v, angular_mom)
        kin_int[idx_molec] = kin_rni[idx_molec] - kin_rot[idx_molec]

        # end loop over molecules

    return {'tot': kin_tot,
            'tra': kin_tra,
            'rni': kin_rni,
            'rot': kin_rot,
            'int': kin_int}


def group_kinetic_energy(kin_molec, nmolecs, molec_group=None):
    r"""
    Sums up the partitioned kinetic energy for a 
    given group or the entire system.

    Parameters
    ----------
    kin_molec : list of dict
        Partitioned kinetic energies per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.
    
    Returns
    -------
    kin : dict
        Dictionary of partitioned kinetic energy for the group.
    """
    #
    kin = {'tot': 0, 'tra': 0, 'rni': 0, 'rot': 0, 'int': 0}
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in kin.keys():
            kin[key] += kin_molec[key][idx_molec]

    return kin


def group_ndof(ndof_molec, nmolecs, molec_group=None):
    r"""
    Sums up the partitioned degrees of freedom for a 
    given group or the entire system.

    Parameters
    ----------
    ndof_molec : list of dict
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.
    
    Returns
    -------
    ndof : dict
        Dictionary of partitioned degrees of freedom for the group.
    """
    #
    ndof = {'tot': 0, 'tra': 0, 'rni': 0, 'rot': 0, 'int': 0}
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in ndof.keys():
            ndof[key] += ndof_molec[idx_molec][key]

    return ndof


def calc_temperatures(kin_molec, ndof_molec, nmolecs, molec_group=None):
    r"""
    Calculates the partitioned temperature for a 
    given group or the entire system.

    Parameters
    ----------
    kin_molec : list of dict
        Partitioned kinetic energies per molecule.
    ndof_molec : list of dict
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    molec_group : iterable
        Indeces of the group to be summed up. None defaults to all molecules
        in the system. Default: None.
    
    Returns
    -------
    temp : dict
        Dictionary of partitioned temperatures for the group.
    """

    kin = group_kinetic_energy(kin_molec, nmolecs, molec_group)
    ndof = group_ndof(ndof_molec, nmolecs, molec_group)

    temp = {}
    for key in kin:
        temp[key] = temperature(kin[key], ndof[key])

    return temp


def test_mb_dist(kin_molec, ndof_molec, nmolecs,
                 temp, alpha, dict_keys, group=None,
                 verbosity=0):
    r"""
    Tests if the partitioned kinetic energy trajectory of a group (or, 
    if group is None, of the entire system) are separately Maxwell-Boltzmann 
    distributed.

    Parameters
    ----------
    kin_molec : list of list of dict
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : list of dict
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    temp : float
        Target temperature of the simulation.
    alpha : float
        Confidence for Maxwell-Boltzmann test.
    dict_keys : list
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    group : iterable
        Indeces of the group to be tested. None defaults to all molecules
        in the system. Default: None.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.
    
    Returns
    -------
    result : int
        Number of partitions for which the Maxwell-Boltzmann distribution
        hypthesis was rejected.
    """
    # save the partitioned kinetic energy trajectories
    group_kin = {key: [] for key in dict_keys}
    ndof = group_ndof(ndof_molec, nmolecs, group)
    # loop over frames
    for k in kin_molec:
        frame_kin = group_kinetic_energy(k, nmolecs, group)
        for key in dict_keys:
            group_kin[key].append(frame_kin[key])

    result = 0
    # test tot, tra, rni, rot, int
    if verbosity > 1:
        print('Testing whether\n'
              '* total,\n'
              '* translational,\n'
              '* rotational & internal,\n'
              '* rotational, and\n'
              '* internal\n'
              'kinetic energies are Maxwell-Boltzmann distributed.')
    elif verbosity > 0:
        print('Testing whether kinetic energies are Maxwell-Boltzmann distributed.')

    for key in dict_keys:
        if check_mb_ensemble(group_kin[key], temp, ndof[key],
                             alpha, verbosity > 2, True, key):
            if verbosity > 1:
                print('* {}: passed'.format(key))
        else:
            result += 1
            if verbosity > 1:
                print('* {}: failed'.format(key))

    if verbosity > 0:
        if result == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    return result


def test_temp_diff(kin_molec, ndof_molec, nmolecs,
                   dtemp, dict_keys, group=None,
                   verbosity=0):
    r"""
    Tests if the partitioned temperatures (averaged over a trajectory)
    of a group (or, if group is None, of the entire system) are within a 
    range `dtemp` of the total temperature.

    Parameters
    ----------
    kin_molec : list of list of dict
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : list of dict
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    dtemp : float
        Target temperature of the simulation.
    dict_keys : list
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    group : iterable
        Indeces of the group to be tested. None defaults to all molecules
        in the system. Default: None.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.
    
    Returns
    -------
    result : int
        Number of partitions for which the temperature ratio to the total
        temperature was larger than dtemp.
    """
    # save the partitioned temperature trajectories
    group_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group)
        for key in dict_keys:
            group_temp[key].append(frame_temp[key])
    # average temperature
    for key in dict_keys:
        group_temp[key] = np.mean(group_temp[key])

    result = 0
    if verbosity > 0:
        print('Testing whether temperatures')
        print('  ' + str(dict_keys[1:]))
        print('are within {:.1f}% of temperature'.format(dtemp*100))
        print('  ' + dict_keys[0] + ' (' + str(group_temp[dict_keys[0]]) + ')')

    temp0 = group_temp[dict_keys[0]]
    for key in dict_keys[1:]:
        if (1 - dtemp) * temp0 <= group_temp[key] <= (1 + dtemp) * temp0:
            if verbosity > 1:
                print('* {} ({:f}): passed'.format(key, group_temp[key]))
        else:
            result += 1
            if verbosity > 1:
                print('* {} ({:f}): failed'.format(key, group_temp[key]))

    if verbosity > 0:
        if result == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    return result


def test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                          group1, group2,
                          dtemp, dict_keys, verbosity=0):
    r"""
    Tests if the partitioned temperatures (averaged over a trajectory) 
    of two groups are (individually) within a range `dtemp` of each others.

    Parameters
    ----------
    kin_molec : list of list of dict
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec : list of dict
        Partitioned degrees of freedom per molecule.
    nmolecs : int
        Total number of molecules in the system.
    group1 : iterable
        Indeces of the first group to be compared.
    group2 : iterable
        Indeces of the second group to be compared.
    dtemp : float
        Target temperature of the simulation.
    dict_keys : list
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    verbosity : int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 0.
        
    Returns
    -------
    result : int
        Number of partitions for which the temperature ratio to the total
        temperature was larger than dtemp.
    """
    # save the partitioned temperature trajectories (group1)
    group1_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group1)
        for key in dict_keys:
            group1_temp[key].append(frame_temp[key])
    # average temperature
    for key in dict_keys:
        group1_temp[key] = np.mean(group1_temp[key])

    # save the partitioned temperature trajectories (group2)
    group2_temp = {key: [] for key in dict_keys}
    # loop over frames
    for k in kin_molec:
        frame_temp = calc_temperatures(k, ndof_molec, nmolecs, group2)
        for key in dict_keys:
            group2_temp[key].append(frame_temp[key])
    # average temperature
    for key in dict_keys:
        group2_temp[key] = np.mean(group2_temp[key])

    result = 0
    if verbosity > 0:
        print('Testing whether temperatures of both groups are within {:.1f}%'.
              format(dtemp*100))

    for key in dict_keys:
        if (1. - dtemp) <= group1_temp[key]/group2_temp[key] <= (1. + dtemp):
            if verbosity > 1:
                print('* {}: passed'.format(key))
        else:
            result += 1
            if verbosity > 1:
                print('* {}: failed'.format(key))

    if verbosity > 0:
        if result == 0:
            print('-> Passed')
        else:
            print('-> Failed')

    return result
