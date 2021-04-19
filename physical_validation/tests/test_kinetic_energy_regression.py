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
import os
from contextlib import redirect_stdout
from io import StringIO
from typing import Dict, List, Optional, Tuple

import pytest

import physical_validation as pv

from .test_system_database import database


def run_kinetic_energy_distribution_check(
    simulation_data: pv.data.SimulationData,
    image_filename_strict: str,
    image_filename_non_strict: str,
) -> Dict:
    r"""

    Parameters
    ----------
    simulation_data
    image_filename_strict
    image_filename_non_strict

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
        bs_seed=23,
    )
    print("\n## Validating kinetic energy distribution (non-strict)")
    result["non-strict"] = pv.kinetic_energy.distribution(
        simulation_data,
        verbosity=2,
        strict=False,
        filename=image_filename_non_strict,
        bs_seed=23,
    )
    return result


def kinetic_energy_distribution_check(
    image_filename_strict: str, image_filename_non_strict: str
) -> Dict:
    system_name = "Water900"
    print("### Regression test of Kinetic energy distribution check using numpy arrays")
    print("### System: " + system_name)

    system = database.system(system_name)
    print("## Creating result object")
    observables = pv.data.ObservableData(
        kinetic_energy=system.observable_as_array("NVT-low", "kinetic_energy"),
        potential_energy=system.observable_as_array("NVT-low", "potential_energy"),
        total_energy=system.observable_as_array("NVT-low", "total_energy"),
        volume=system.observable_as_array("NVT-low", "volume"),
    )
    simulation_data = pv.data.SimulationData(
        units=system.units,
        ensemble=system.ensemble("NVT-low"),
        system=system.system_data,
        observables=observables,
    )

    return run_kinetic_energy_distribution_check(
        simulation_data,
        image_filename_strict=image_filename_strict,
        image_filename_non_strict=image_filename_non_strict,
    )


def test_kinetic_energy_distribution_regression(
    data_regression, file_regression, image_regression
) -> None:
    test_output = StringIO()
    test_image_strict = "test_plot_strict.png"
    test_image_non_strict = "test_plot_non_strict.png"

    def remove_image(test_image):
        # Remove image if it exists
        try:
            os.remove(test_image)
        except OSError:
            pass

    remove_image(test_image_strict)
    remove_image(test_image_non_strict)

    # Redirect stdout into string which we can test
    with redirect_stdout(test_output):
        result = kinetic_energy_distribution_check(
            test_image_strict, test_image_non_strict
        )

    # Regression fixture can only check dicts of strings
    # Use 6 digits to keep reproducibility reasonable
    result["strict"] = "{:.6f}".format(result["strict"])
    result["non-strict"] = (
        "(" + ",".join(["{:.6f}".format(v) for v in result["non-strict"]]) + ")"
    )
    data_regression.check(result)
    # Test printed output
    file_regression.check(contents=test_output.getvalue())
    # Test printed pictures
    with open(test_image_strict, "rb") as image:
        image_regression.check(
            image_data=image.read(),
            diff_threshold=1.0,
            basename="test_kinetic_energy_distribution_regression_strict",
        )
    with open(test_image_non_strict, "rb") as image:
        image_regression.check(
            image_data=image.read(),
            diff_threshold=1.0,
            basename="test_kinetic_energy_distribution_regression_non_strict",
        )

    # Clean up
    remove_image(test_image_strict)
    remove_image(test_image_non_strict)
