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
This module contains low-level functionality of the
`physical_validation.kinetic_energy` module. The functions in this module
should generally not be called directly. Please use the high-level
functions from `physical_validation.kinetic energy`.
"""
import warnings
from multiprocessing.pool import ThreadPool
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import scipy.stats as stats

from ..util import trajectory
from . import plot


def is_close(
    value1: float,
    value2: float,
    relative_tolerance: float = 1e-09,
    absolute_tolerance: float = 1e-09,
) -> bool:
    r"""
    Whether two float values are close, as defined by a given relative and
    absolute tolerance.
    """
    return abs(value1 - value2) <= max(
        relative_tolerance * max(abs(value1), abs(value2)), absolute_tolerance
    )


def check_distribution(
    kin: np.ndarray,
    temp: float,
    ndof: float,
    kb: float,
    verbosity: int,
    screen: bool,
    filename: Optional[str],
    ene_unit: Optional[str],
    temp_unit: Optional[str],
    data_is_uncorrelated: bool,
) -> float:
    r"""
    Checks if a kinetic energy trajectory is Maxwell-Boltzmann distributed.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using
       the high-level version based on the SimulationData object. See
       physical_validation.kinetic_energy.check_mb_ensemble for more
       information and full documentation.

    Parameters
    ----------
    kin
        Kinetic energy snapshots of the system.
    temp
        Target temperature of the system. Used to construct the
        Maxwell-Boltzmann distribution.
    ndof
        Number of degrees of freedom in the system. Used to construct the
        Maxwell-Boltzmann distribution.
    kb
        Boltzmann constant :math:`k_B`.
    verbosity
        0: Silent.
        1: Print minimal information.
        2: Print result details.
        3: Print additional information.
    screen
        Plot distributions on screen.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    ene_unit
        Energy unit - used for output only.
    temp_unit
        Temperature unit - used for output only.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result
        The p value of the test.

    See Also
    --------
    physical_validation.kinetic_energy.distribution : High-level version
    """

    # Discard burn-in period and time-correlated frames
    kin = trajectory.prepare(
        kin,
        verbosity=verbosity,
        name="Kinetic energy",
        skip_preparation=data_is_uncorrelated,
    )
    kt = kb * temp

    if ndof <= 0:
        warnings.warn("Zero degrees of freedom!")
        p = float("NaN")
    else:
        d, p = stats.kstest(kin, "gamma", (ndof / 2, 0, kt))

    # ====================== #
    # Plot to screen or file #
    # ====================== #
    do_plot = screen or filename is not None
    if do_plot:
        ana_dist = stats.gamma(ndof / 2, scale=kt)
        ana_kin = np.linspace(ana_dist.ppf(0.0001), ana_dist.ppf(0.9999), 200)
        ana_hist = ana_dist.pdf(ana_kin)

        tunit = ""
        if temp_unit is not None:
            tunit = temp_unit

        data = [
            {
                "y": kin,
                "hist": True,
                "args": dict(label="Trajectory", density=True, alpha=0.5),
            }
        ]
        if ndof > 0:
            data.append(
                {
                    "x": ana_kin,
                    "y": ana_hist,
                    "args": dict(label="Analytical T=" + str(temp) + tunit, lw=5),
                }
            )

        unit = ""
        if ene_unit is not None:
            unit = " [" + ene_unit + "]"

        plot.plot(
            data,
            legend="lower left",
            title="Kinetic energy distribution",
            xlabel="Kinetic energy" + unit,
            ylabel="Probability [%]",
            sci_x=True,
            percent=True,
            filename=filename,
            screen=screen,
        )

    if verbosity > 0:
        if verbosity > 1:
            message = (
                "Kinetic energy distribution check (strict)\n"
                "Kolmogorov-Smirnov test result: p = {:g}\n"
                "Null hypothesis: Kinetic energy is Maxwell-Boltzmann distributed".format(
                    p
                )
            )
        else:
            message = "p = {:g}".format(p)
        print(message)

    return p


def check_mean_std(
    kin: np.ndarray,
    temp: float,
    ndof: float,
    kb: float,
    verbosity: int,
    bs_repetitions: int,
    bootstrap_seed: Optional[int],
    screen: bool,
    filename: Optional[str],
    ene_unit: Optional[str],
    temp_unit: Optional[str],
    data_is_uncorrelated: bool,
) -> Tuple[float, float]:
    r"""
    Calculates the mean and standard deviation of a trajectory (+ bootstrap
    error estimates), and compares them to the theoretically expected values.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using
       the high-level version based on the SimulationData object. See
       physical_validation.kinetic_energy.check_mb_ensemble for more
       information and full documentation.

    Parameters
    ----------
    kin
        Kinetic energy snapshots of the system.
    temp
        Target temperature of the system. Used to construct the
        Maxwell-Boltzmann distribution.
    ndof
        Number of degrees of freedom in the system. Used to construct the
        Maxwell-Boltzmann distribution.
    kb
        Boltzmann constant :math:`k_B`.
    verbosity
        0: Silent.
        1: Print minimal information.
        2: Print result details.
        3: Print additional information.
    bs_repetitions
        Number of bootstrap samples used for error estimate.
    bootstrap_seed
        Sets the random number seed for bootstrapping.
        If set, bootstrapping will be reproducible.
        If `None`, bootstrapping is non-reproducible.
    screen
        Plot distributions on screen.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    ene_unit
        Energy unit - used for output only.
    temp_unit
        Temperature unit - used for output only.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result
        Distance of the estimated T(mu) and T(sigma) from the expected
        temperature, measured in standard deviations of the estimates.

    See Also
    --------
    physical_validation.kinetic_energy.distribution : High-level version
    """

    # Discard burn-in period and time-correlated frames
    kin = trajectory.prepare(
        kin,
        verbosity=verbosity,
        name="Kinetic energy",
        skip_preparation=data_is_uncorrelated,
    )

    if ndof <= 0:
        warnings.warn("Zero degrees of freedom!")

    # ========================== #
    # Compute mu and sig of data #
    # ========================== #
    kt = temp * kb
    loc = 0
    ana_shape = ndof / 2
    ana_scale = kt
    ana_dist = stats.gamma(ana_shape, loc=loc, scale=ana_scale)

    if ndof > 0:
        temp_mu = 2 * np.mean(kin) / (ndof * kb)
        temp_sig = np.sqrt(2 / ndof) * np.std(kin) / kb
    else:
        temp_mu = 0
        temp_sig = 0

    # ======================== #
    # Bootstrap error estimate #
    # ======================== #
    if bootstrap_seed is not None:
        np.random.seed(bootstrap_seed)
    mu = []
    sig = []
    for k in trajectory.bootstrap(kin, bs_repetitions):
        mu.append(np.mean(k))
        sig.append(np.std(k))
    std_mu = np.std(mu)
    std_sig = np.std(sig)
    if ndof > 0:
        std_temp_mu = 2 * std_mu / (ndof * kb)
        std_temp_sig = np.sqrt(2 / ndof) * std_sig / kb
    else:
        std_temp_mu = 0
        std_temp_sig = 0

    # ====================== #
    # Plot to screen or file #
    # ====================== #
    do_plot = screen or filename is not None
    if do_plot:
        ana_kin = np.linspace(ana_dist.ppf(0.0001), ana_dist.ppf(0.9999), 200)
        ana_hist = ana_dist.pdf(ana_kin)

        tunit = ""
        if temp_unit is not None:
            tunit = temp_unit

        data = [
            {
                "y": kin,
                "hist": True,
                "args": dict(label="Trajectory", density=True, alpha=0.5),
            }
        ]
        if ndof > 0:
            data.append(
                {
                    "x": ana_kin,
                    "y": ana_hist,
                    "args": dict(label="Analytical T=" + str(temp) + tunit, lw=5),
                }
            )

        unit = ""
        if ene_unit is not None:
            unit = " [" + ene_unit + "]"

        plot.plot(
            data,
            legend="best",
            title="Kinetic energy distribution",
            xlabel="Kinetic energy" + unit,
            ylabel="Probability [%]",
            sci_x=True,
            percent=True,
            filename=filename,
            screen=screen,
        )

    # ================ #
    # Output to screen #
    # ================ #
    if verbosity > 0:
        eunit = ""
        if ene_unit is not None:
            eunit = " " + ene_unit
        tunit = ""
        if temp_unit is not None:
            tunit = " " + temp_unit
        if verbosity > 1:
            message = (
                "Kinetic energy distribution check (non-strict)\n"
                "Analytical distribution (T={2:.2f}{0:s}):\n"
                " * mu: {3:.2f}{1:s}\n"
                " * sigma: {4:.2f}{1:s}\n"
                "Trajectory:\n"
                " * mu: {5:.2f} +- {7:.2f}{1:s}\n"
                "   T(mu) = {9:.2f} +- {11:.2f}{0:s}\n"
                " * sigma: {6:.2f} +- {8:.2f}{1:s}\n"
                "   T(sigma) = {10:.2f} +- {12:.2f}{0:s}".format(
                    tunit,
                    eunit,
                    temp,
                    ana_dist.mean(),
                    ana_dist.std(),
                    np.mean(kin),
                    np.std(kin),
                    std_mu,
                    std_sig,
                    temp_mu,
                    temp_sig,
                    std_temp_mu,
                    std_temp_sig,
                )
            )
        else:
            message = (
                "T(mu) = {1:.2f} +- {3:.2f}{0:s}\n"
                "T(sigma) = {2:.2f} +- {4:.2f}{0:s}".format(
                    tunit, temp_mu, temp_sig, std_temp_mu, std_temp_sig
                )
            )
        print(message)

    # ============= #
    # Return values #
    # ============= #
    nan = float("NaN")
    if ndof > 0:
        r1 = np.abs(temp - temp_mu) / std_temp_mu
        r2 = np.abs(temp - temp_sig) / std_temp_sig
    else:
        r1 = nan
        r2 = nan
    return r1, r2


def check_equipartition(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    molec_idx: np.ndarray,
    molec_nbonds: np.ndarray,
    natoms: int,
    nmolecs: int,
    temp: float,
    kb: float,
    strict: bool,
    ndof_reduction_tra: float,
    ndof_reduction_rot: float,
    molec_groups: Optional[List[np.ndarray]],
    random_divisions: int,
    random_groups: int,
    random_division_seed: Optional[int],
    ndof_molec: Optional[List[Dict[str, float]]],
    kin_molec: Optional[List[List[Dict[str, float]]]],
    verbosity: int,
    screen: bool,
    filename: Optional[str],
    ene_unit: Optional[str],
    temp_unit: Optional[str],
    bootstrap_seed: Optional[int],
    data_is_uncorrelated: bool,
) -> Tuple[
    Union[List[float], List[Tuple[float, float]]],
    List[Dict[str, float]],
    List[List[Dict[str, float]]],
]:
    r"""
    Checks the equipartition of a simulation trajectory.

    .. warning: This is a low-level function. Additionally to being less
       user-friendly, there is a higher probability of erroneous and / or
       badly documented behavior due to unexpected inputs. Consider using
       the high-level version based on the SimulationData object. See
       physical_validation.kinetic_energy.check_equipartition for more
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
    natoms
        Number of atoms in the system
    nmolecs
        Number of molecules in the system
    temp
        Target temperature of the simulation.
    kb
        Boltzmann constant :math:`k_B`.
    strict
        If True, check full kinetic energy distribution via K-S test.
        Otherwise, check mean and width of kinetic energy distribution.
    ndof_reduction_tra
        Number of center-of-mass translational degrees of freedom to
        remove.
    ndof_reduction_rot
        Number of center-of-mass rotational degrees of freedom to remove.
    molec_groups : List[array-like] (ngroups x ?), optional
        List of 1d arrays containing molecule indeces defining groups.
        Useful to pre-define groups of molecules (e.g. solute / solvent,
        liquid mixture species, ...). If None, no pre-defined molecule
        groups will be tested.

        *Note:* If an empty 1d array is found as last element in the list, the remaining
        molecules are collected in this array. This allows, for example, to only
        specify the solute, and indicate the solvent by giving an empty array.
    random_divisions
        Number of random division tests attempted.
    random_groups
        Number of groups the system is randomly divided in.
    random_division_seed
        Seed making the random divisions reproducible.
        If `None`, random divisions not reproducible
    ndof_molec
        Pass in the degrees of freedom per molecule. Slightly increases speed of repeated
        analysis of the same simulation run.
    kin_molec
        Pass in the kinetic energy per molecule. Greatly increases speed of repeated
        analysis of the same simulation run.
    verbosity
        Verbosity level, where 0 is quiet and 3 very chatty.
    screen
        Plot distributions on screen. Default: False.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    ene_unit
        Energy unit - used for output only.
    temp_unit
        Temperature unit - used for output only.
    bootstrap_seed
        Sets the random number seed for bootstrapping (if strict=False).
        If set, bootstrapping will be reproducible.
        If `None`, bootstrapping is non-reproducible.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result : List[float] or List[Tuple[float]]
        If `strict=True`: The p value for every test.
        If `strict=False`: Distance of the estimated T(mu) and T(sigma) from
        the expected temperature, measured in standard deviations of the
        respective estimate, for every test.
    ndof_molec : List[dict]
        List of the degrees of freedom per molecule. Can be saved to increase speed of
        repeated analysis of the same simulation run.
    kin_molec : List[List[dict]]
        List of the kinetic energy per molecule per frame. Can be saved to increase speed
        of repeated analysis of the same simulation run.

    See Also
    --------
    physical_validation.kinetic_energy.check_equipartition : High-level version

    """

    dict_keys = [
        "total",
        "translational",
        "rotational and internal",
        "rotational",
        "internal",
    ]
    result = []

    # ========================== #
    # Compute molecular energies #
    # ========================== #
    # for each molecule, calculate total / translational / rotational & internal /
    #   rotational / internal degrees of freedom
    #   returns: List[dict] of floats (shape: nmolecs x 5 x 1)
    if ndof_molec is None:
        ndof_molec = calc_ndof(
            natoms,
            nmolecs,
            molec_idx,
            molec_nbonds,
            ndof_reduction_tra,
            ndof_reduction_rot,
        )

    # for each frame, calculate total / translational / rotational & internal /
    #   rotational / internal kinetic energy for each molecule
    if kin_molec is None:
        with ThreadPool() as p:
            kin_molec = p.starmap(
                calc_molec_kinetic_energy,
                [
                    (r, v, masses, molec_idx, natoms, nmolecs)
                    for r, v in zip(positions, velocities)
                ],
            )

    # =============================== #
    # Check system-wide equipartition #
    # =============================== #
    result.extend(
        test_group(
            kin_molec=kin_molec,
            ndof_molec=ndof_molec,
            nmolecs=nmolecs,
            temp=temp,
            kb=kb,
            dict_keys=dict_keys,
            strict=strict,
            group=None,
            verbosity=verbosity,
            screen=screen,
            filename=filename,
            ene_unit=ene_unit,
            temp_unit=temp_unit,
            bootstrap_seed=bootstrap_seed,
            data_is_uncorrelated=data_is_uncorrelated,
        )
    )

    # ==================================== #
    # Check equipartition in random groups #
    # ==================================== #
    if random_division_seed is not None:
        np.random.seed(random_division_seed)
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
                print("Testing randomly divided group {:d}".format(rg))
            if verbosity > 3:
                print(group)
            result.extend(
                test_group(
                    kin_molec=kin_molec,
                    ndof_molec=ndof_molec,
                    nmolecs=nmolecs,
                    temp=temp,
                    kb=kb,
                    dict_keys=dict_keys,
                    strict=strict,
                    group=group,
                    verbosity=verbosity,
                    screen=screen,
                    filename=filename,
                    ene_unit=ene_unit,
                    temp_unit=temp_unit,
                    bootstrap_seed=bootstrap_seed,
                    data_is_uncorrelated=data_is_uncorrelated,
                )
            )

    # ======================================== #
    # Check equipartition in predefined groups #
    # ======================================== #
    # use predefined group division?
    # if no groups, return
    if not molec_groups:
        return result, ndof_molec, kin_molec
    # is last group empty?
    last_empty = not molec_groups[-1]
    # get rid of empty groups
    molec_groups = [g for g in molec_groups if g]
    # if no groups now (all were empty), return now
    if not molec_groups:
        return result, ndof_molec, kin_molec

    if last_empty:
        # last group is [] -> insert remaining molecules
        combined = []
        for group in molec_groups:
            combined.extend(group)
        molec_groups.append(np.array([m for m in range(nmolecs) if m not in combined]))

    for mg, group in enumerate(molec_groups):
        if verbosity > 0:
            print("Testing predifined divided group {:d}".format(mg))
        if verbosity > 3:
            print(group)
        result.extend(
            test_group(
                kin_molec=kin_molec,
                ndof_molec=ndof_molec,
                nmolecs=nmolecs,
                temp=temp,
                kb=kb,
                dict_keys=dict_keys,
                strict=strict,
                group=group,
                verbosity=verbosity,
                screen=screen,
                filename=filename,
                ene_unit=ene_unit,
                temp_unit=temp_unit,
                bootstrap_seed=bootstrap_seed,
                data_is_uncorrelated=data_is_uncorrelated,
            )
        )

    return result, ndof_molec, kin_molec


def calc_ndof(
    natoms: int,
    nmolecs: int,
    molec_idx: np.ndarray,
    molec_nbonds: np.ndarray,
    ndof_reduction_tra: float,
    ndof_reduction_rot: float,
) -> List[Dict[str, float]]:
    r"""
    Calculates the total / translational / rotational & internal /
    rotational / internal degrees of freedom per molecule.

    Parameters
    ----------
    natoms
        Total number of atoms in the system
    nmolecs
        Total number of molecules in the system
    molec_idx
        Index of first atom for every molecule
    molec_nbonds
        Number of bonds for every molecule
    ndof_reduction_tra
        Number of center-of-mass translational degrees of freedom tobremove.
    ndof_reduction_rot
        Number of center-of-mass rotational degrees of freedom to remove.

    Returns
    -------
    ndof_molec
        List of dictionaries containing the degrees of freedom for each molecule
        Keys: ['total', 'translational', 'rotational and internal', 'rotational', 'internal']
    """
    # check whether there are monoatomic molecules:
    nmono = (np.array(molec_idx[1:]) - np.array(molec_idx[:-1]) == 1).sum()

    # ndof to be deducted per molecule
    # ndof reduction due to COM motion constraining
    ndof_com_tra_pm = ndof_reduction_tra / nmolecs
    ndof_com_rot_pm = ndof_reduction_rot / (nmolecs - nmono)

    ndof_molec = []
    # add last idx to molec_idx to ease looping
    molec_idx = np.append(molec_idx, [natoms])
    # loop over molecules
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(
        zip(molec_idx[:-1], molec_idx[1:])
    ):
        natoms = idx_atm_end - idx_atm_init
        nbonds = molec_nbonds[idx_molec]
        ndof_tot = 3 * natoms - nbonds - ndof_com_tra_pm - ndof_com_rot_pm
        ndof_tra = 3 - ndof_com_tra_pm
        ndof_rni = ndof_tot - ndof_tra
        ndof_rot = 3 - ndof_com_rot_pm
        ndof_int = ndof_tot - ndof_tra - ndof_rot
        if is_close(ndof_int, 0, absolute_tolerance=1e-09):
            ndof_int = 0
        if natoms == 1:
            ndof_tot = 3 - ndof_com_tra_pm
            ndof_tra = 3 - ndof_com_tra_pm
            ndof_rni = 0
            ndof_rot = 0
            ndof_int = 0
        ndof_molec.append(
            {
                "total": ndof_tot,
                "translational": ndof_tra,
                "rotational and internal": ndof_rni,
                "rotational": ndof_rot,
                "internal": ndof_int,
            }
        )

    return ndof_molec


def calc_molec_kinetic_energy(
    pos: np.ndarray,
    vel: np.ndarray,
    masses: np.ndarray,
    molec_idx: np.ndarray,
    natoms: int,
    nmolecs: int,
) -> Dict[str, np.ndarray]:
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
    natoms
        Total number of atoms in the system
    nmolecs
        Total number of molecules in the system

    Returns
    -------
    kin
        Dictionary of lists containing the kinetic energies for each molecule
        Keys: ['total', 'translational', 'rotational and internal', 'rotational', 'internal']
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
    for idx_molec, (idx_atm_init, idx_atm_end) in enumerate(
        zip(molec_idx[:-1], molec_idx[1:])
    ):
        # if monoatomic molecule
        if idx_atm_end == idx_atm_init + 1:
            v = vel[idx_atm_init]
            m = masses[idx_atm_init]
            kin_tot[idx_molec] = 0.5 * m * np.dot(v, v)
            kin_tra[idx_molec] = kin_tot[idx_molec]
            kin_rni[idx_molec] = 0
            kin_rot[idx_molec] = 0
            kin_int[idx_molec] = 0
            continue
        # compute center of mass position, velocity and total mass
        com_r = np.zeros(3)
        com_v = np.zeros(3)
        com_m = 0
        # loop over atoms in molecule
        for r, v, m in zip(
            pos[idx_atm_init:idx_atm_end],
            vel[idx_atm_init:idx_atm_end],
            masses[idx_atm_init:idx_atm_end],
        ):
            com_r += m * r
            com_v += m * v
            com_m += m

            # total kinetic energy is straightforward
            kin_tot[idx_molec] += 0.5 * m * np.dot(v, v)

        com_r /= com_m
        com_v /= com_m

        # translational kinetic energy
        kin_tra[idx_molec] = 0.5 * com_m * np.dot(com_v, com_v)
        # combined rotational and internal kinetic energy
        kin_rni[idx_molec] = kin_tot[idx_molec] - kin_tra[idx_molec]

        # compute tensor of inertia and angular momentum
        inertia = np.zeros((3, 3))
        angular_mom = np.zeros(3)
        # loop over atoms in molecule
        for r, v, m in zip(
            pos[idx_atm_init:idx_atm_end],
            vel[idx_atm_init:idx_atm_end],
            masses[idx_atm_init:idx_atm_end],
        ):
            # relative positions and velocities
            rr = r - com_r
            rv = v - com_v
            rr2 = np.dot(rr, rr)
            # inertia tensor:
            #   (i,i) = m*(r*r - r(i)*r(i))
            #   (i,j) = m*r(i)*r(j) (i != j)
            atm_inertia = -m * np.tensordot(rr, rr, axes=0)
            for i in range(3):
                atm_inertia[i][i] += m * rr2
            inertia += atm_inertia
            # angular momentum: r x p
            angular_mom += m * np.cross(rr, rv)

        # angular velocity of the molecule: inertia^{-1} * angular_mom
        angular_v = np.dot(np.linalg.inv(inertia), angular_mom)

        kin_rot[idx_molec] = 0.5 * np.dot(angular_v, angular_mom)
        kin_int[idx_molec] = kin_rni[idx_molec] - kin_rot[idx_molec]

        # end loop over molecules

    return {
        "total": kin_tot,
        "translational": kin_tra,
        "rotational and internal": kin_rni,
        "rotational": kin_rot,
        "internal": kin_int,
    }


def group_kinetic_energy(
    kin_molec: Dict[str, np.ndarray], nmolecs: int, molec_group=Optional[Iterable[int]]
) -> Dict[str, float]:
    r"""
    Sums up the partitioned kinetic energy for a
    given group or the entire system.

    Parameters
    ----------
    kin_molec
        Partitioned kinetic energies per molecule.
    nmolecs
        Total number of molecules in the system.
    molec_group
        Indices of the group to be summed up. `None` defaults to all molecules
        in the system.

    Returns
    -------
    kin : dict
        Dictionary of partitioned kinetic energy for the group.
    """
    #
    kin = {
        "total": 0,
        "translational": 0,
        "rotational and internal": 0,
        "rotational": 0,
        "internal": 0,
    }
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in kin.keys():
            kin[key] += kin_molec[key][idx_molec]

    return kin


def group_ndof(
    ndof_molec: List[Dict[str, float]],
    nmolecs: int,
    molec_group=Optional[Iterable[int]],
) -> Dict[str, float]:
    r"""
    Sums up the partitioned degrees of freedom for a
    given group or the entire system.

    Parameters
    ----------
    ndof_molec
        Partitioned degrees of freedom per molecule.
    nmolecs
        Total number of molecules in the system.
    molec_group
        Indeces of the group to be summed up. None defaults to all molecules
        in the system.

    Returns
    -------
    ndof : dict
        Dictionary of partitioned degrees of freedom for the group.
    """
    #
    ndof = {
        "total": 0,
        "translational": 0,
        "rotational and internal": 0,
        "rotational": 0,
        "internal": 0,
    }
    #
    if molec_group is None:
        molec_group = range(nmolecs)
    # loop over molecules
    for idx_molec in molec_group:
        for key in ndof.keys():
            ndof[key] += ndof_molec[idx_molec][key]

    return ndof


def test_group(
    kin_molec: List[Dict[str, np.ndarray]],
    ndof_molec: List[Dict[str, float]],
    nmolecs: int,
    temp: float,
    kb: float,
    dict_keys: List[str],
    strict: bool,
    group: Optional[Iterable[int]],
    verbosity: int,
    screen: bool,
    filename: Optional[str],
    ene_unit: Optional[str],
    temp_unit: Optional[str],
    bootstrap_seed: Optional[int],
    data_is_uncorrelated: bool,
) -> Union[List[float], List[Tuple[float, float]]]:
    r"""
    Tests if the partitioned kinetic energy trajectory of a group (or,
    if group is None, of the entire system) are separately Maxwell-Boltzmann
    distributed.

    Parameters
    ----------
    kin_molec
        Partitioned kinetic energies per molecule for every frame.
    ndof_molec
        Partitioned degrees of freedom per molecule.
    nmolecs
        Total number of molecules in the system.
    temp
        Target temperature of the simulation.
    kb
        Boltzmann constant :math:`k_B`.
    dict_keys
        List of dictionary keys representing the partitions of the degrees
        of freedom.
    strict
        If True, check full kinetic energy distribution via K-S test.
        Otherwise, check mean and width of kinetic energy distribution.
    group
        Indices of the group to be tested. `None` defaults to all molecules
        in the system.
    verbosity
        Verbosity level, where 0 is quiet and 3 very chatty.
    screen
        Plot distributions on screen.
    filename
        Plot distributions to `filename`. If `None`, no plotting.
    ene_unit
        Energy unit - used for output only.
    temp_unit
        Temperature unit - used for output only.
    bootstrap_seed
        Sets the random number seed for bootstrapping (if strict=False).
        If set, bootstrapping will be reproducible.
        If `None`, bootstrapping is non-reproducible.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result
        p value for every partition (strict) or tuple of distance of
        estimated T(mu) and T(sigma) for every partition (non-strict)
    """
    # save the partitioned kinetic energy trajectories
    group_kin = {key: [] for key in dict_keys}
    ndof = group_ndof(ndof_molec, nmolecs, group)
    # loop over frames
    for k in kin_molec:
        frame_kin = group_kinetic_energy(k, nmolecs, group)
        for key in dict_keys:
            group_kin[key].append(frame_kin[key])

    result = []
    # test total, translational, rotational and internal, rotational, int
    if verbosity > 0:
        message = "Equipartition: Testing group-wise kinetic energies ("
        if not strict:
            message += "non-"
        message += "strict)"
        print(message)

    for key in dict_keys:
        if verbosity > 1:
            print("* " + key + ":")
        if filename is None:
            fn = None
        else:
            fn = key.replace(" ", "_") + "_" + filename
        if ndof[key] < 1e-9:
            print("No {:s} degrees of freedom in this group".format(key))
            continue
        if strict:
            res = check_distribution(
                kin=group_kin[key],
                temp=temp,
                ndof=ndof[key],
                kb=kb,
                verbosity=int(verbosity > 1),
                screen=screen,
                filename=fn,
                ene_unit=ene_unit,
                temp_unit=temp_unit,
                data_is_uncorrelated=data_is_uncorrelated,
            )
        else:
            res = check_mean_std(
                kin=group_kin[key],
                temp=temp,
                ndof=ndof[key],
                kb=kb,
                verbosity=int(verbosity > 1),
                screen=screen,
                filename=fn,
                ene_unit=ene_unit,
                temp_unit=temp_unit,
                bs_repetitions=200,
                bootstrap_seed=bootstrap_seed,
                data_is_uncorrelated=data_is_uncorrelated,
            )
        result.append(res)

    return result
