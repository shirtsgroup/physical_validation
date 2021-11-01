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
The `kinetic_energy` module is part of the physical_validation package, and
consists of checks of the kinetic energy distribution and its
equipartition.
"""
from typing import List, Optional, Tuple, Union

import numpy as np

from .data import ObservableData, SimulationData
from .util import kinetic_energy as util_kin


def distribution(
    data: SimulationData,
    strict: bool = False,
    verbosity: int = 2,
    screen: bool = False,
    filename: Optional[str] = None,
    bs_repetitions: int = 200,
    bootstrap_seed: Optional[int] = None,
    data_is_uncorrelated: bool = False,
) -> Union[float, Tuple[float, float]]:
    r"""Checks the distribution of a kinetic energy trajectory.

    Parameters
    ----------
    data
        Simulation data object
    strict
        If True, check full kinetic energy distribution via K-S test.
        Otherwise, check mean and width of kinetic energy distribution.
        Default: False
    verbosity
        Verbosity level, where 0 is quiet and 3 shows full details. Default: 2.
    screen
        Plot distributions on screen. Default: False.
    filename
        Plot distributions to `filename`. Default: None, no plotting to file.
    bs_repetitions
        Number of bootstrap samples used for error estimate (if strict=False).
        Default: 200.
    bootstrap_seed
        Sets the random number seed for bootstrapping (if strict=False).
        If set, bootstrapping will be reproducible.
        Default: `None`, bootstrapping is non-reproducible.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result
        If `strict=True`: The p value of the test.
        If `strict=False`: Distance of the estimated T(mu) and T(sigma) from
        the expected temperature, measured in standard deviations of the respective
        estimate.

    Notes
    -----
    Non-strict test
        If `strict = False` (the default), this function will estimate the mean and
        the standard deviation of the data. Analytically, a gamma distribution with
        shape :math:`k = N / 2` (with :math:`N` the number of degrees of freedom)
        and scale :math:`\theta = k_B T` (with :math:`T` the target temperature)
        is expected. The mean and the standard deviation of a gamma distribution
        are given by :math:`\mu = k\theta` and :math:`\sigma = \sqrt k \theta`.

        The standard error of the mean and standard deviation are estimated via
        bootstrap resampling. The function prints the analytically expected mean
        and variance as well as the fitted values and their error estimates. It
        also prints T(mu) and T(sigma), which are defined as the temperatures to
        which the estimated mean and standard deviation correspond, given the number
        of degrees of freedom :math:`N` in the system:

        .. math::
            T(\mu') = \frac{2 \mu'}{N k_B}

        .. math::
            T(\sigma') = \frac{\sqrt 2 \sigma'}{\sqrt N k_B}

        The return value is a tuple containing the distance of the estimated T(mu) and
        T(sigma) from the expected temperature, measured in standard deviations of the
        respective estimates.

    Strict test
        If `strict = True`, this function tests the hypothesis that a sample
        of kinetic energies is Maxwell-Boltzmann distributed given a specific
        target temperature and the number of degrees of freedom in the system,

        .. math::
            P(K) \sim K^{N/2-1} e^{-\beta K} \, .

        The test is performed using the Kolmogorov-Smirnov test provided by
        scipy.stats.kstest_. It returns the :math:`p`-value, measuring the
        likelihood that a sample at least as extreme as the one given is
        originating from the expected distribution.

        .. _scipy.stats.kstest: https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.kstest.html

        .. note:: The Kolmogorov-Smirnov test is known to have two weaknesses.

           #. The test is more sensitive towards deviations around the center
              of the distribution than at its tails. We deem this to be acceptable
              for most MD applications, but be wary if yours is sensible to the
              kinetic distribution tails.
           #. The test is not valid if its parameters are guessed from the data
              set. Using the target temperature of the MD simulation as an input
              is therefore perfectly valid, but using the average temperature
              over the trajectory as an input to the test can potentially
              invalidate it.

    """
    data.raise_if_units_are_none(
        test_name="kinetic_energy.distribution",
        argument_name="data",
    )
    data.raise_if_ensemble_is_invalid(
        test_name="kinetic_energy.distribution",
        argument_name="data",
        check_pressure=False,
        check_mu=False,
    )
    data.raise_if_system_data_is_invalid(
        test_name="kinetic_energy.distribution",
        argument_name="data",
        check_full_system_data_only=True,
    )
    data.raise_if_observable_data_is_invalid(
        required_observables=["kinetic_energy"],
        test_name="kinetic_energy.distribution",
        argument_name="data",
    )

    ndof = (
        data.system.natoms * 3
        - data.system.nconstraints
        - data.system.ndof_reduction_tra
        - data.system.ndof_reduction_rot
    )

    if strict:
        return util_kin.check_distribution(
            kin=data.observables.kinetic_energy,
            temp=data.ensemble.temperature,
            ndof=ndof,
            kb=data.units.kb,
            verbosity=verbosity,
            screen=screen,
            filename=filename,
            ene_unit=data.units.energy_str,
            temp_unit=data.units.temperature_str,
            data_is_uncorrelated=data_is_uncorrelated,
        )
    else:
        return util_kin.check_mean_std(
            kin=data.observables.kinetic_energy,
            temp=data.ensemble.temperature,
            ndof=ndof,
            kb=data.units.kb,
            verbosity=verbosity,
            bs_repetitions=bs_repetitions,
            bootstrap_seed=bootstrap_seed,
            screen=screen,
            filename=filename,
            ene_unit=data.units.energy_str,
            temp_unit=data.units.temperature_str,
            data_is_uncorrelated=data_is_uncorrelated,
        )


def equipartition(
    data: SimulationData,
    strict: bool = False,
    molec_groups: Optional[List[np.ndarray]] = None,
    random_divisions: int = 0,
    random_groups: int = 0,
    random_division_seed: Optional[int] = None,
    verbosity: int = 2,
    screen: bool = False,
    filename: Optional[str] = None,
    bootstrap_seed: Optional[int] = None,
    data_is_uncorrelated: bool = False,
) -> Union[List[float], List[Tuple[float, float]]]:
    r"""Checks the equipartition of a simulation trajectory.

    Parameters
    ----------
    data
        Simulation data object
    strict
        If True, check full kinetic energy distribution via K-S test.
        Otherwise, check mean and width of kinetic energy distribution.
        Default: False
    molec_groups
        List of 1d arrays containing molecule indices defining groups. Useful to pre-define
        groups of molecules (e.g. solute / solvent, liquid mixture species, ...). If None,
        no pre-defined molecule groups will be tested. Default: None.

        *Note:* If an empty 1d array is found as last element in the list, the remaining
        molecules are collected in this array. This allows, for example, to only
        specify the solute, and indicate the solvent by giving an empty array.
    random_divisions
        Number of random division tests attempted. Default: 0 (random division tests off).
    random_groups
        Number of groups the system is randomly divided in. Default: 2.
    random_division_seed
        Seed making the random divisions reproducible.
        Default: None, random divisions not reproducible
    verbosity
        Verbosity level, where 0 is quiet and 3 very chatty. Default: 2.
    screen
        Plot distributions on screen. Default: False.
    filename
        Plot distributions to `filename`. Default: None, no plotting to file
    bootstrap_seed
        Sets the random number seed for bootstrapping (if strict=False).
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping is non-reproducible.
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    result
        If `strict=True`: The p value for every tests.
        If `strict=False`: Distance of the estimated T(mu) and T(sigma) from
        the expected temperature, measured in standard deviations of the
        respective estimate, for every test.

    Notes
    -----
    This function compares the kinetic energy between groups of degrees of
    freedom. Theoretically, the kinetic energy is expected (via the
    equipartition theorem) to be equally distributed over all degrees of
    freedom. In practice, deviations of temperature between groups of
    degrees of freedom up to several degrees K are routinely observed.
    Larger deviations can, however, hint to misbehaving simulations, such
    as, e.g., frozen degrees of freedom, lack of energy exchange between
    degrees of freedom, and transfer of heat from faster to slower degrees
    of freedom.

    Splitting of degrees of freedom is done both on a sub-molecular and on
    a molecular level. On a sub-molecular level, the degrees of freedom of
    a molecule can be partitioned into rigid-body contributions
    (translation of the center-of-mass, rotation around the
    center-of-mass) and intra-molecular contributions. On a molecular
    level, the single molecules of the system can be divided in groups,
    either by function (solute / solvent, different species of liquid
    mixtures, ...) or randomly.

    `check_equipartition()` partitions the kinetic energy of the entire
    system and, optionally, of predefined or randomly separated groups. It
    then computes either the mean and the standard deviation of each partition
    and compares them to the theoretically expected value (`strict=True`, the
    default), or it performs a Kolmogorov-Smirnov test of the distribution.
    See physical_validation.kinetic_energy.distribution for more detail on the
    checks.

    """
    if data.observables is None:
        # The equipartition test doesn't need observable input, but uses the
        # data structure to store information. Create it if it doesn't exist.
        data.observables = ObservableData()

    data.raise_if_units_are_none(
        test_name="kinetic_energy.equipartition",
        argument_name="data",
    )
    data.raise_if_ensemble_is_invalid(
        test_name="kinetic_energy.equipartition",
        argument_name="data",
        check_pressure=False,
        check_mu=False,
    )
    data.raise_if_system_data_is_invalid(
        test_name="kinetic_energy.equipartition",
        argument_name="data",
        check_full_system_data_only=False,
    )
    data.raise_if_trajectory_data_is_invalid(
        test_name="kinetic_energy.equipartition",
        argument_name="data",
    )

    (
        result,
        data.system.ndof_per_molecule,
        data.observables.kinetic_energy_per_molecule,
    ) = util_kin.check_equipartition(
        positions=data.trajectory["position"],
        velocities=data.trajectory["velocity"],
        masses=data.system.mass,
        molec_idx=data.system.molecule_idx,
        molec_nbonds=data.system.nconstraints_per_molecule,
        natoms=data.system.natoms,
        nmolecs=len(data.system.molecule_idx),
        temp=data.ensemble.temperature,
        kb=data.units.kb,
        strict=strict,
        ndof_reduction_tra=data.system.ndof_reduction_tra,
        ndof_reduction_rot=data.system.ndof_reduction_rot,
        molec_groups=molec_groups,
        random_divisions=random_divisions,
        random_groups=random_groups,
        random_division_seed=random_division_seed,
        ndof_molec=data.system.ndof_per_molecule,
        kin_molec=data.observables.kinetic_energy_per_molecule,
        verbosity=verbosity,
        screen=screen,
        filename=filename,
        ene_unit=data.units.energy_str,
        temp_unit=data.units.temperature_str,
        bootstrap_seed=bootstrap_seed,
        data_is_uncorrelated=data_is_uncorrelated,
    )

    return result
