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
This software allows users to perform statistical test to determine if a
given molecular simulation is consistent with the thermodynamic ensemble
it is performed in.

Users should cite the JCTC paper: Shirts, M. R. "Simple Quantitative
Tests to Validate Sampling from Thermodynamic Ensembles",
J. Chem. Theory Comput., 2013, 9 (2), pp 909-926,
https://dx.doi.org/10.1021/ct300688p
"""
from typing import Dict, List, Optional

import numpy as np

from .data import SimulationData
from .util import ensemble
from .util import error as pv_error


def beta_warning() -> None:
    print(
        "###############################################################################"
    )
    print(
        "# WARNING: Support for muVT ensemble is an experimental feature under current #"
    )
    print(
        "#          development. You can help us to improve it by reporting errors     #"
    )
    print(
        "#          at https://github.com/shirtsgroup/physical_validation              #"
    )
    print(
        "#          Thank you!                                                         #"
    )
    print(
        "###############################################################################"
    )


def check(
    data_sim_one: SimulationData,
    data_sim_two: SimulationData,
    total_energy: bool = False,
    bootstrap_error: bool = False,
    bootstrap_repetitions: int = 200,
    bootstrap_seed: Optional[int] = None,
    screen: bool = False,
    filename: Optional[str] = None,
    verbosity: int = 1,
    data_is_uncorrelated: bool = False,
) -> List[float]:
    r"""
    Check the ensemble. The correct check is inferred from the
    simulation data given.

    Parameters
    ----------
    data_sim_one
        Simulation data object of first simulation
    data_sim_two
        Simulation data object of second simulation differing in
        its state point from the first
    total_energy
        Whether to use the total energy for the calculation
        Default: False, use potential energy only
    bootstrap_error
        Calculate the standard error via bootstrap resampling
        Default: False
    bootstrap_repetitions
        Number of bootstrap repetitions drawn
        Default: 200
    bootstrap_seed
        Sets the random number seed for bootstrapping.
        If set, bootstrapping will be reproducible.
        Default: None, bootstrapping is non-reproducible.
    screen
        Plot distributions on screen. Default: False.
    filename
        Plot distributions to `filename`.
        Default: `None`, no plotting to file.
    verbosity
        Level of verbosity, from 0 (quiet) to 3 (very verbose).
        Default: 1
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.

    Returns
    -------
    quantiles
        The number of quantiles the computed result is off the analytical one.

    """
    data_sim_one.raise_if_units_are_none(
        test_name="ensemble.check",
        argument_name="data_sim_one",
    )
    data_sim_two.raise_if_units_are_none(
        test_name="ensemble.check",
        argument_name="data_sim_two",
    )
    if not SimulationData.compatible(data_sim_one, data_sim_two):
        raise pv_error.InputError(
            ["data_sim_one", "data_sim_two"], "Simulation data not compatible."
        )
    data_sim_one.raise_if_ensemble_is_invalid(
        test_name="ensemble.check",
        argument_name="data_sim_one",
        check_pressure=True,
        check_mu=True,
    )
    data_sim_two.raise_if_ensemble_is_invalid(
        test_name="ensemble.check",
        argument_name="data_sim_two",
        check_pressure=True,
        check_mu=True,
    )

    if data_sim_one.ensemble.ensemble != data_sim_two.ensemble.ensemble:
        raise pv_error.InputError(
            ["data_sim_one", "data_sim_two"],
            "The two simulations were sampling different ensembles. "
            "The simulations are expected to differ in state point "
            "(e.g. target temperature, target pressure), but not "
            "in their sampled ensemble (e.g. NVT, NPT).",
        )

    sampled_ensemble = data_sim_one.ensemble.ensemble

    if not (
        sampled_ensemble == "NVT"
        or sampled_ensemble == "NPT"
        or sampled_ensemble == "muVT"
    ):
        raise NotImplementedError(
            "Test of ensemble " + sampled_ensemble + " is not implemented (yet).",
        )

    labels = {
        "E": "Total Energy",
        "U": "Potential Energy",
        "H": "Enthalpy",
        "V": "Volume",
        "N": "Number of Species",
        r"E - \sum \mu N": r"$E - \sum \mu N$",
        r"U - \sum \mu N": r"$U - \sum \mu N$",
    }

    energy_observable_name = "total_energy" if total_energy else "potential_energy"
    energy_observable_abbreviation = "E" if total_energy else "U"

    energy_units = data_sim_one.units.energy_str

    quantiles = None

    num_bins_for_analysis = 40
    tail_cutoff = 0.001  # 0.1%

    if sampled_ensemble == "NVT":
        data_sim_one.raise_if_observable_data_is_invalid(
            required_observables=[energy_observable_name],
            test_name="ensemble.check",
            argument_name="data_sim_one",
        )
        data_sim_two.raise_if_observable_data_is_invalid(
            required_observables=[energy_observable_name],
            test_name="ensemble.check",
            argument_name="data_sim_two",
        )
        quantiles = ensemble.check_1d(
            traj1=data_sim_one.observables[energy_observable_name],
            traj2=data_sim_two.observables[energy_observable_name],
            param1=data_sim_one.ensemble.temperature,
            param2=data_sim_two.ensemble.temperature,
            kb=data_sim_one.units.kb,
            quantity=energy_observable_abbreviation,
            dtemp=True,
            dpress=False,
            dmu=False,
            temp=None,
            pvconvert=None,
            nbins=num_bins_for_analysis,
            cutoff=tail_cutoff,
            bootstrap_seed=bootstrap_seed,
            bootstrap_error=bootstrap_error,
            bootstrap_repetitions=bootstrap_repetitions,
            verbosity=verbosity,
            filename=filename,
            screen=screen,
            xlabel=labels[energy_observable_abbreviation],
            xunit=energy_units,
            data_is_uncorrelated=data_is_uncorrelated,
        )

    elif sampled_ensemble == "NPT":
        temperatures = np.array(
            [data_sim_one.ensemble.temperature, data_sim_two.ensemble.temperature]
        )
        pressures = np.array(
            [data_sim_one.ensemble.pressure, data_sim_two.ensemble.pressure]
        )
        equal_temps = temperatures[0] == temperatures[1]
        equal_press = pressures[0] == pressures[1]

        required_observables = ["volume"]
        if equal_press:
            required_observables.append(energy_observable_name)

        data_sim_one.raise_if_observable_data_is_invalid(
            required_observables=required_observables,
            test_name="ensemble.check",
            argument_name="data_sim_one",
        )
        data_sim_two.raise_if_observable_data_is_invalid(
            required_observables=required_observables,
            test_name="ensemble.check",
            argument_name="data_sim_two",
        )

        volume_units = data_sim_one.units.volume_str

        # Calculate conversion from p*V to energy units
        #
        # GROMACS standard units are
        #   energy: kJ/mol
        #   volume: nm^3
        #   pressure: bar
        #   => pV-term: bar * nm^3 == 1e-25 kJ == 6.022140857e-2 kJ/mol
        #   => pvconvert = 6.022140857e-2
        # UnitData stores conversion factors relative to GROMACS units
        #   energy: energy_conversion * kJ/mol
        #   volume: volume_conversion * nm^3
        #   pressure: pressure_conversion * bar
        #   => pV-term: [p]*[V] == pressure_conversion * volume_conversion bar * nm^3
        # Units were checked earlier, so we can use either simulation data structure
        pvconvert = 6.022140857e-2
        pvconvert *= (
            data_sim_one.units.pressure_conversion
            * data_sim_one.units.volume_conversion
        )
        pvconvert /= data_sim_one.units.energy_conversion

        if equal_press and not equal_temps:
            if energy_observable_abbreviation == "U":
                energy_observable_abbreviation = "H"
            quantiles = ensemble.check_1d(
                traj1=data_sim_one.observables[energy_observable_name]
                + pvconvert * pressures[0] * data_sim_one.observables.volume,
                traj2=data_sim_two.observables[energy_observable_name]
                + pvconvert * pressures[1] * data_sim_two.observables.volume,
                param1=temperatures[0],
                param2=temperatures[1],
                kb=data_sim_one.units.kb,
                quantity=energy_observable_abbreviation,
                dtemp=True,
                dpress=False,
                dmu=False,
                temp=None,
                pvconvert=None,
                nbins=num_bins_for_analysis,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                xlabel=labels[energy_observable_abbreviation],
                xunit=energy_units,
                data_is_uncorrelated=data_is_uncorrelated,
            )
        elif equal_temps and not equal_press:
            quantiles = ensemble.check_1d(
                traj1=data_sim_one.observables.volume,
                traj2=data_sim_two.observables.volume,
                param1=pressures[0],
                param2=pressures[1],
                kb=data_sim_one.units.kb,
                quantity="V",
                dtemp=False,
                dpress=True,
                dmu=False,
                nbins=num_bins_for_analysis,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                temp=temperatures[0],
                pvconvert=pvconvert,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                xlabel=labels["V"],
                xunit=volume_units,
                data_is_uncorrelated=data_is_uncorrelated,
            )
        else:
            param1 = np.array([temperatures[0], pressures[0]])
            param2 = np.array([temperatures[1], pressures[1]])
            quantiles = ensemble.check_2d(
                traj1=np.array(
                    [
                        data_sim_one.observables[energy_observable_name],
                        data_sim_one.observables.volume,
                    ]
                ),
                traj2=np.array(
                    [
                        data_sim_two.observables[energy_observable_name],
                        data_sim_two.observables.volume,
                    ]
                ),
                param1=param1,
                param2=param2,
                kb=data_sim_one.units.kb,
                pvconvert=pvconvert,
                quantity=[energy_observable_abbreviation, "V"],
                dtempdpress=True,
                dtempdmu=False,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                data_is_uncorrelated=data_is_uncorrelated,
            )

    elif sampled_ensemble == "muVT":
        beta_warning()

        equal_temps = (
            data_sim_one.ensemble.temperature == data_sim_two.ensemble.temperature
        )
        num_differing_chemical_potentials = np.count_nonzero(
            data_sim_one.ensemble.mu != data_sim_two.ensemble.mu
        )
        num_identical_chemical_potentials = np.count_nonzero(
            data_sim_one.ensemble.mu == data_sim_two.ensemble.mu
        )

        if num_differing_chemical_potentials > 1:
            raise NotImplementedError(
                "For muVT, ensemble checking is only implemented for chemical potentials "
                "differing for at most one species."
            )

        if not equal_temps and num_differing_chemical_potentials == 0:
            # quantity = U - S_i mu_i*N_i
            # parameter = T
            energy_observable_abbreviation += r" - \sum \mu N"
            quantiles = ensemble.check_1d(
                traj1=data_sim_one.observables[energy_observable_name]
                - ensemble.chemical_potential_energy(
                    data_sim_one.ensemble.mu, data_sim_one.observables.number_of_species
                ),
                traj2=data_sim_two.observables[energy_observable_name]
                - ensemble.chemical_potential_energy(
                    data_sim_two.ensemble.mu, data_sim_two.observables.number_of_species
                ),
                param1=data_sim_one.ensemble.temperature,
                param2=data_sim_two.ensemble.temperature,
                kb=data_sim_one.units.kb,
                quantity=energy_observable_abbreviation,
                dtemp=True,
                dpress=False,
                dmu=False,
                temp=None,
                pvconvert=None,
                nbins=num_bins_for_analysis,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                xlabel=labels[energy_observable_abbreviation],
                xunit=energy_units,
                data_is_uncorrelated=data_is_uncorrelated,
            )

        if equal_temps and num_differing_chemical_potentials == 1:
            # quantity = N_i
            # parameter = mu_i
            idx = np.flatnonzero(data_sim_one.ensemble.mu != data_sim_two.ensemble.mu)[
                0
            ]
            quantiles = ensemble.check_1d(
                traj1=data_sim_one.observables.number_of_species[:, idx],
                traj2=data_sim_two.observables.number_of_species[:, idx],
                param1=data_sim_one.ensemble.mu[idx],
                param2=data_sim_two.ensemble.mu[idx],
                kb=data_sim_one.units.kb,
                quantity="N",
                dtemp=False,
                dpress=False,
                dmu=True,
                nbins=num_bins_for_analysis,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                temp=data_sim_one.ensemble.temperature,
                pvconvert=None,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                xlabel=labels["N"],
                xunit=None,
                data_is_uncorrelated=data_is_uncorrelated,
            )

        if not equal_temps and num_differing_chemical_potentials == 1:
            # quantity = U - S_i!=j mu_i*N_i | mu_j*N_j
            # parameter = beta, mu_j
            energy_observable_one = data_sim_one.observables[energy_observable_name]
            energy_observable_two = data_sim_two.observables[energy_observable_name]
            if num_identical_chemical_potentials > 0:
                idx_identical = np.flatnonzero(
                    data_sim_one.ensemble.mu == data_sim_two.ensemble.mu
                )
                energy_observable_one -= ensemble.chemical_potential_energy(
                    data_sim_one.ensemble.mu[idx_identical],
                    data_sim_one.observables.number_of_species[:, idx_identical],
                )
                energy_observable_two -= ensemble.chemical_potential_energy(
                    data_sim_two.ensemble.mu[idx_identical],
                    data_sim_two.observables.number_of_species[:, idx_identical],
                )
                energy_observable_abbreviation += r" - \sum mu*N"

            idx_differing = np.flatnonzero(
                data_sim_one.ensemble.mu != data_sim_two.ensemble.mu
            )[0]

            param1 = np.array(
                [
                    data_sim_one.ensemble.temperature,
                    data_sim_one.ensemble.mu[idx_differing],
                ]
            )
            param2 = np.array(
                [
                    data_sim_two.ensemble.temperature,
                    data_sim_two.ensemble.mu[idx_differing],
                ]
            )
            quantiles = ensemble.check_2d(
                traj1=np.array(
                    [
                        energy_observable_one,
                        data_sim_one.observables.number_of_species[:, idx_differing],
                    ]
                ),
                traj2=np.array(
                    [
                        energy_observable_two,
                        data_sim_two.observables.number_of_species[:, idx_differing],
                    ]
                ),
                param1=param1,
                param2=param2,
                kb=data_sim_one.units.kb,
                pvconvert=None,
                quantity=[energy_observable_abbreviation, "N"],
                dtempdpress=False,
                dtempdmu=True,
                cutoff=tail_cutoff,
                bootstrap_seed=bootstrap_seed,
                bootstrap_error=bootstrap_error,
                bootstrap_repetitions=bootstrap_repetitions,
                verbosity=verbosity,
                filename=filename,
                screen=screen,
                data_is_uncorrelated=data_is_uncorrelated,
            )

    return quantiles


def estimate_interval(
    data: SimulationData,
    verbosity: int = 1,
    total_energy: bool = False,
    data_is_uncorrelated: bool = False,
) -> Dict[str, float]:
    r"""
    In order to perform an ensemble check, two simulations at distinct state
    point are needed. Choosing two state points too far apart will result
    in poor or zero overlap between the distributions, leading to very noisy
    results (due to sample errors in the tails) or a breakdown of the method,
    respectively. Choosing two state points very close to each others, on the
    other hand, makes it difficult to distinguish the slope from statistical
    error in the samples.

    This function implements a rule of thumb based on the standard deviations
    of distributions. It takes a single simulation and suggests appropriate
    intervals for a second simulation to be used for ensemble checking.

    Parameters
    ----------
    data
        The performed simulation.
    verbosity
        If 0, no output is printed on screen. If 1, estimated intervals are
        printed. If larger, additional information during calculation are
        printed.
        Default: 1
    total_energy
        Use total energy instead of potential energy only.
        Default: False
    data_is_uncorrelated
        Whether the provided data is uncorrelated. If this option
        is set, the equilibration, decorrelation and tail pruning
        of the trajectory is skipped. This can speed up the analysis,
        but note that if the provided data is correlated, the results
        of the physical validation checks might be invalid.
        Default: False

    Returns
    -------
    intervals
        If `data` was performed under NVT conditions, `intervals` contains only
        one entry:

            * `'dT'`, containing the suggested temperature interval.

        If `data` was performed under NPT conditions, `intervals` contains three
        entries:

            * `'dT'`: Suggested temperature interval at constant pressure
            * `'dP'`: Suggested pressure interval at constant temperature
            * `'dTdP'`: Suggested combined temperature and pressure interval

    """
    data.raise_if_units_are_none(
        test_name="ensemble.estimate_interval",
        argument_name="data",
    )
    data.raise_if_ensemble_is_invalid(
        test_name="ensemble.estimate_interval",
        argument_name="data",
        check_pressure=True,
        check_mu=True,
    )
    if not (
        data.ensemble.ensemble == "NVT"
        or data.ensemble.ensemble == "NPT"
        or data.ensemble.ensemble == "muVT"
    ):
        raise NotImplementedError(
            "estimate_interval() not implemented for ensemble " + data.ensemble.ensemble
        )

    energy_observable_name = "total_energy" if total_energy else "potential_energy"
    required_observables = [energy_observable_name]
    if data.ensemble.ensemble == "NPT":
        required_observables.append("volume")
    if data.ensemble.ensemble == "muVT":
        required_observables.append("number_of_species")
    data.raise_if_observable_data_is_invalid(
        required_observables=required_observables,
        test_name="ensemble.estimate_interval",
        argument_name="data",
    )

    tail_cutoff = 0.001  # 0.1%

    if data.ensemble.ensemble == "NVT":
        result = ensemble.estimate_interval(
            ens_string="NVT",
            ens_temp=data.ensemble.temperature,
            energy=data.observables[energy_observable_name],
            kb=data.units.kb,
            ens_press=None,
            volume=None,
            pvconvert=None,
            ens_mu=None,
            species_number=None,
            verbosity=verbosity,
            cutoff=tail_cutoff,
            tunit=data.units.temperature_str,
            punit="",
            munit="",
            data_is_uncorrelated=data_is_uncorrelated,
        )
    elif data.ensemble.ensemble == "NPT":
        pvconvert = 6.022140857e-2
        pvconvert *= data.units.pressure_conversion * data.units.volume_conversion
        pvconvert /= data.units.energy_conversion
        result = ensemble.estimate_interval(
            ens_string="NPT",
            ens_temp=data.ensemble.temperature,
            energy=data.observables[energy_observable_name],
            kb=data.units.kb,
            ens_press=data.ensemble.pressure,
            volume=data.observables.volume,
            pvconvert=pvconvert,
            ens_mu=None,
            species_number=None,
            verbosity=verbosity,
            cutoff=tail_cutoff,
            tunit=data.units.temperature_str,
            punit=data.units.pressure_str,
            munit="",
            data_is_uncorrelated=data_is_uncorrelated,
        )
    else:
        beta_warning()

        result = ensemble.estimate_interval(
            ens_string="muVT",
            ens_temp=data.ensemble.temperature,
            energy=data.observables[energy_observable_name],
            kb=data.units.kb,
            ens_press=None,
            volume=None,
            pvconvert=None,
            ens_mu=data.ensemble.mu,
            species_number=data.observables.number_of_species,
            verbosity=verbosity,
            cutoff=tail_cutoff,
            tunit=data.units.temperature_str,
            punit="",
            munit=data.units.energy_str,
            data_is_uncorrelated=data_is_uncorrelated,
        )

    return result
