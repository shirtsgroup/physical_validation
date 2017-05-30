from __future__ import division
from __future__ import print_function

import numpy as np

import physicalvalidation.util.kineticenergy as util_eq
from physicalvalidation.deprecated.check_input import equipartition


def check_equipartition(positions, velocities, masses,
                        molec_idx, molec_nbonds, molec_groups=None,
                        dtemp=0.1, temp=None, alpha=0.05,
                        stop_com_tra=False, stop_com_rot=False,
                        random_divisions=0, random_groups=2, verbosity=2):
    """
    Checks the equipartition of a simulation trajectory.
    
    This function compares the kinetic energy between groups of degrees of
    freedom. Theoretically, the kinetic energy is expected (via the equipartition 
    theorem) to be equally distributed over all degrees of freedom. In practice, 
    deviations of temperature between groups of degrees of freedom up to several 
    degrees K are routinely observed. Larger deviations can, however, hint to 
    misbehaving simulations, such as, e.g., frozen degrees of freedom, lack of
    energy exchange between degrees of freedom, and transfer of heat from faster
    to slower oscillating degrees of freedom.
    
    Splitting of degrees of freedom is done both on a sub-molecular and on a 
    molecular level. On a sub-molecular level, the degrees of freedom of a molecule
    can be partitioned into rigid-body contributions (translation of the center-of-mass,
    rotation around the center-of-mass) and intra-molecular contributions. On a 
    molecular level, the single molecules of the system can be divided in groups, 
    either by function (solute / solvent, different species of liquid mixtures, ...)
    or randomly.
    
    `check_equipartition` compares the partitioned temperatures of the entire system
    and, optionally, of predefined or randomly separated groups. 
    
    Note: In theory, the kinetic energy of the subgroups are expected to be individually 
    Maxwell-Boltzmann distributed. As this is seldomly holding in practice (see above),
    `check_equipartition` is by default checking only for abnormal deviations in average
    temperatures. The more strict Maxwell-Boltzmann testing can be invoked by giving the
    target temperature `temp` as an input.
    
    TODO: Should it actually be "no external forces" instead of "was the center-of-mass
          removed during simulation?"?
    
    :param positions: nd-array (nframes x natoms x 3)
        3d array containing the positions of all atoms for all frames
    :param velocities: nd-array (nframes x natoms x 3)
        3d array containing the velocities of all atoms for all frames
    :param masses: nd-array (natoms x 1)
        1d array containing the masses of all atoms
    :param molec_idx: nd-array (nmolecs x 1)
        Index of first atom for every molecule 
    :param molec_nbonds: nd-array (nmolecs x 1)
        Number of bonds for every molecule
    :param molec_groups: list of nd-array (ngroups x ?) 
        List of 1d arrays containing molecule indeces defining groups. Useful to pre-define
        groups of molecules (e.g. solute / solvent, liquid mixture species, ...). If None,
        no pre-defined molecule groups will be tested. Default: None.
    :param dtemp: float
        Fraction of temperature deviation tolerated between groups. Default: 0.05 (5%).
    :param temp: float
        Target temperature of the simulation. If None, the kinetic energies will not be
        tested for Maxwell-Boltzmann distribution, but only compared amongst each others.
        Default: None.
    :param alpha: float
        Confidence for Maxwell-Boltzmann test. Default: 0.05 (5%).
    :param stop_com_tra: bool
        Was the center-of-mass translation removed during the simulation? Default: False.
    :param stop_com_rot: bool
        Was the center-of-mass translation removed during the simulation? Default: False.
    :param random_divisions: int
        Number of random division tests attempted. Default: 0 (random division tests off).
    :param random_groups: int
        Number of groups the system is randomly divided in. Default: 2.
    :param verbosity: int
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 2.
    :return: int
        Number of equipartition violations. Tune up verbosity for details.
    """

    dict_keys = ['tot', 'tra', 'rni', 'rot', 'int']

    # do input validation, calculate number of atoms / molecules / frames
    # TODO: Check this. Accept array-like, take a more relaxed approach to input validation.
    natoms, nmolecs, nframes = equipartition(positions, velocities, masses,
                                             molec_idx, molec_nbonds)

    # for each molecule, calculate total / translational / rotational & internal /
    #   rotational / internal degrees of freedom
    #   returns: list of dict of floats (shape: nmolecs x 5 x 1)
    ndof_molec = util_eq.calc_ndof(natoms, nmolecs, molec_idx, molec_nbonds,
                                   stop_com_tra, stop_com_rot)

    # for each frame, calculate total / translational / rotational & internal /
    #   rotational / internal kinetic energy for each molecule
    kin_molec = []
    for r, v in zip(positions, velocities):
        kin_molec.append(util_eq.calc_kinetic_energy(r, v, masses,
                                                     molec_idx, natoms, nmolecs))

    result = 0

    # test system-wide tot, tra, rni, rot, int
    if temp is not None:
        # check for Maxwell-Boltzmann distribution of
        # partitioned kinetic energy trajectories
        result += util_eq.test_mb_dist(kin_molec, ndof_molec, nmolecs,
                                       temp, alpha, dict_keys, verbosity)
    else:
        # compare partitioned temperatures to total temperature
        result += util_eq.test_temp_diff(kin_molec, ndof_molec, nmolecs,
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
                result += util_eq.test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                               alpha, dict_keys, group, verbosity)
            else:
                result += util_eq.test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                                 dtemp, dict_keys, group, verbosity)
        # test groups against each others
        for rg1, group1 in enumerate(groups):
            for group2 in groups[rg1+1:]:
                util_eq.test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                              group1, group2,
                                              dtemp, dict_keys, verbosity)

    # use predefined division
    if molec_groups is not None:
        for mg, group in enumerate(molec_groups):
            if verbosity > 0:
                print('Testing predifined divided group {:d}'.format(mg))
            if verbosity > 3:
                print(group)
            if temp is not None:
                result += util_eq.test_mb_dist(kin_molec, ndof_molec, nmolecs, temp,
                                               alpha, dict_keys, group, verbosity)
            else:
                result += util_eq.test_temp_diff(kin_molec, ndof_molec, nmolecs,
                                                 dtemp, dict_keys, group, verbosity)
        # test groups against each others
        for rg1, group1 in enumerate(molec_groups):
            for group2 in molec_groups[rg1+1:]:
                util_eq.test_temp_diff_groups(kin_molec, ndof_molec, nmolecs,
                                              group1, group2,
                                              dtemp, dict_keys, verbosity)

    return result
