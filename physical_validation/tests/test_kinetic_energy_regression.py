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
Regression tests running the kinetic energy checks and making sure
that the output is unchanged between code versions.

Output changes might be intended. If you are sure a change is wanted,
you can update the reference data using
    $ pytest test_kinetic_energy_regression.py --force-regen
from within the `physical_validation/tests` folder and commit the
changes. Note that due to the different regression test fixtures
used (data, text and image), you might need to run the above command
repeatedly to update all reference data.

The regression tests require `pytest-regressions`
(https://pytest-regressions.readthedocs.io).
"""
import glob
import os
from contextlib import redirect_stdout
from io import StringIO
from typing import Dict

import pytest

import physical_validation as pv

from .test_system_database import database


def run_kinetic_energy_distribution_check(
    simulation_data: pv.data.SimulationData,
    image_filename_strict: str,
    image_filename_non_strict: str,
) -> Dict:
    r"""
    Runs the kinetic energy distribution check in its strict and non-strict
    variant on the provided simulation data object.

    Parameters
    ----------
    simulation_data
        Simulation data object
    image_filename_strict
        File name to plot the results of the strict test
    image_filename_non_strict
        File name to plot the results of the non-strict test

    Returns
    -------
    Dict containing the return values of the strict and non-strict tests
    """

    result = {}
    print("\n## Validating kinetic energy distribution (strict)")
    result["strict"] = pv.kinetic_energy.distribution(
        simulation_data,
        verbosity=2,
        strict=True,
        filename=image_filename_strict,
        bootstrap_seed=23,
    )
    print("\n## Validating kinetic energy distribution (non-strict)")
    result["non-strict"] = pv.kinetic_energy.distribution(
        simulation_data,
        verbosity=2,
        strict=False,
        filename=image_filename_non_strict,
        bootstrap_seed=23,
    )

    # Regression fixture can only check dicts of strings
    # Use 6 digits to keep reproducibility reasonable
    result["strict"] = "{:.6f}".format(result["strict"])
    result["non-strict"] = (
        "(" + ",".join(["{:.6f}".format(v) for v in result["non-strict"]]) + ")"
    )
    return result


def run_kinetic_energy_equipartition_check(
    simulation_data: pv.data.SimulationData,
    image_filename_strict: str,
    image_filename_non_strict: str,
) -> Dict:
    r"""

    Parameters
    ----------
    simulation_data
        Simulation data object
    image_filename_strict
        File name to plot the results of the strict test
    image_filename_non_strict
        File name to plot the results of the non-strict test

    Returns
    -------
    Dict containing the return values of the strict and non-strict tests
    """

    result = {}
    print("\n## Validating kinetic energy equipartition (strict)")
    # Run once calculating only the rotational, translational and internal
    # degrees of freedom for the entire system, and plotting the results
    pv.kinetic_energy.equipartition(
        simulation_data,
        strict=True,
        molec_groups=None,
        random_divisions=0,
        random_groups=0,
        verbosity=2,
        filename=image_filename_strict,
        bootstrap_seed=7,
        data_is_uncorrelated=True,
    )
    # Run again including molec_groups and random divisions, but no plots
    # to speed up execution
    result["strict"] = pv.kinetic_energy.equipartition(
        simulation_data,
        strict=True,
        molec_groups=[[0], []],
        random_divisions=1,
        random_groups=2,
        random_division_seed=14,
        verbosity=2,
        bootstrap_seed=7,
        data_is_uncorrelated=True,
    )
    print("\n## Validating kinetic energy equipartition (non-strict)")
    # Run once calculating only the rotational, translational and internal
    # degrees of freedom for the entire system, and plotting the results
    pv.kinetic_energy.equipartition(
        simulation_data,
        strict=False,
        molec_groups=None,
        random_divisions=0,
        random_groups=0,
        verbosity=2,
        filename=image_filename_non_strict,
        bootstrap_seed=7,
        data_is_uncorrelated=True,
    )
    # Run again including molec_groups and random divisions, but no plots
    # to speed up execution
    result["non-strict"] = pv.kinetic_energy.equipartition(
        simulation_data,
        strict=False,
        molec_groups=[[0], []],
        random_divisions=1,
        random_groups=2,
        random_division_seed=19,
        verbosity=2,
        bootstrap_seed=7,
        data_is_uncorrelated=True,
    )

    # Regression fixture can only check dicts of strings
    # Use 6 digits to keep reproducibility reasonable
    result["strict"] = (
        "[" + ",".join(["{:.6f}".format(v) for v in result["strict"]]) + "]"
    )
    result["non-strict"] = (
        "["
        + ",".join(
            [
                "(" + ",".join(["{:.6f}".format(v) for v in res]) + ")"
                for res in result["non-strict"]
            ]
        )
        + "]"
    )

    return result


def kinetic_energy_distribution_check(
    image_filename_strict: str, image_filename_non_strict: str
) -> Dict:
    system_name = "Water900"
    print("### Regression test of kinetic energy distribution check using numpy arrays")
    print("### System: " + system_name)

    system = database.system(system_name)
    print("## Creating result object")
    parser = pv.data.FlatfileParser()
    simulation_data = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble("NVT-low"),
        system=system.system_data,
        kinetic_ene_file=system.observable_flat_file("NVT-low", "kinetic_energy"),
        potential_ene_file=system.observable_flat_file("NVT-low", "potential_energy"),
        total_ene_file=system.observable_flat_file("NVT-low", "total_energy"),
        volume_file=system.observable_flat_file("NVT-low", "volume"),
    )

    return run_kinetic_energy_distribution_check(
        simulation_data,
        image_filename_strict=image_filename_strict,
        image_filename_non_strict=image_filename_non_strict,
    )


def kinetic_energy_equipartition_check(
    image_filename_strict: str, image_filename_non_strict: str
) -> Dict:
    system_name = "Octanol2"
    print("### Regression test of kinetic energy equipartition check using flat files")
    print("### System: " + system_name)

    system = database.system(system_name)
    print("## Creating result object")
    parser = pv.data.FlatfileParser()
    simulation_data = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble("GasPhase"),
        system=system.system_data,
        position_file=system.trajectory_flat_file("GasPhase", "position"),
        velocity_file=system.trajectory_flat_file("GasPhase", "velocity"),
    )

    return run_kinetic_energy_equipartition_check(
        simulation_data,
        image_filename_strict=image_filename_strict,
        image_filename_non_strict=image_filename_non_strict,
    )


# We have small test trajectories, so ignore plotting warnings
@pytest.mark.filterwarnings("ignore:For smooth histograms")
@pytest.mark.parametrize("test_type", ["distribution", "equipartition"])
def test_kinetic_energy_regression(
    data_regression,
    file_regression,
    image_regression,
    test_type: str,
) -> None:
    test_output = StringIO()
    test_image_strict = f"{test_type}_plot_strict.png"
    test_image_non_strict = f"{test_type}_plot_non_strict.png"
    test_name = f"test_kinetic_energy_regression_{test_type}"

    def remove_image(test_image):
        if test_type == "distribution":
            # Remove image if it exists
            try:
                os.remove(test_image)
            except OSError:
                pass
        elif test_type == "equipartition":
            # Remove images if they exist
            for file in glob.glob(f"*{test_image}"):
                os.remove(file)

    def check_image_base(test_image):
        with open(test_image, "rb") as image:
            image_regression.check(
                image_data=image.read(),
                diff_threshold=1.0,
                basename=f"{test_name}_{test_image}",
            )

    def check_image(test_image):
        if test_type == "distribution":
            check_image_base(test_image)
        elif test_type == "equipartition":
            for file in glob.glob(f"*{test_image}"):
                check_image_base(file)

    remove_image(test_image_strict)
    remove_image(test_image_non_strict)

    # Redirect stdout into string which we can test
    with redirect_stdout(test_output):
        if test_type == "distribution":
            result = kinetic_energy_distribution_check(
                image_filename_strict=test_image_strict,
                image_filename_non_strict=test_image_non_strict,
            )
        elif test_type == "equipartition":
            result = kinetic_energy_equipartition_check(
                image_filename_strict=test_image_strict,
                image_filename_non_strict=test_image_non_strict,
            )
        else:
            raise NotImplementedError(f"Unknown test type {test_type}")

    # Test returned values
    data_regression.check(result, basename=test_name)
    # Test printed output
    file_regression.check(contents=test_output.getvalue(), basename=test_name)
    # Test printed pictures
    check_image(test_image_strict)
    check_image(test_image_non_strict)

    # Clean up
    remove_image(test_image_strict)
    remove_image(test_image_non_strict)
