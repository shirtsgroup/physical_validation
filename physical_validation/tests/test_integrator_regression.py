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
Regression tests running the integrator convergence checks and making sure
that the output is unchanged between code versions.

Output changes might be intended. If you are sure a change is wanted,
you can update the reference data using
    $ pytest test_integrator_regression.py --force-regen
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
from typing import Dict, List

import physical_validation as pv

from .test_system_database import database


def run_integrator_convergence(
    simulations: List[pv.data.SimulationData], image_filename: str
) -> Dict[str, str]:
    result = pv.integrator.convergence(simulations=simulations, filename=image_filename)

    # Regression fixture can only check dicts of strings
    # Use 6 digits to keep reproducibility reasonable
    return {"result": "{:.6f}".format(result)}


def integrator_convergence(image_filename: str) -> Dict[str, str]:
    system_name = "Argon1000"
    print("### Regression test of integrator convergence check")
    print("### System: " + system_name)

    parser = pv.data.FlatfileParser()
    system = database.system(system_name)
    simulation_data_list = []
    for time_step in system.simulations:
        simulation_data_list.append(
            parser.get_simulation_data(
                units=system.units,
                ensemble=system.ensemble("liquid"),
                system=system.system_data,
                const_of_mot_file=system.observable_flat_file(
                    time_step, "constant of motion"
                ),
                dt=float(time_step),
            )
        )

    return run_integrator_convergence(
        simulations=simulation_data_list, image_filename=image_filename
    )


def test_integrator_convergence(
    data_regression, file_regression, image_regression
) -> None:
    r"""
    TODO: something

    Parameters
    ----------
    data_regression
        Regression test fixture testing python dicts
    file_regression
        Regression test fixture testing text files
    image_regression
        Regression test fixture testing images
    """
    test_output = StringIO()
    test_image = "test_plot.png"

    # Remove image if it exists
    try:
        os.remove(test_image)
    except OSError:
        pass

    # Redirect stdout into string which we can test
    with redirect_stdout(test_output):
        result = integrator_convergence(test_image)

    data_regression.check(result)
    # Test printed output
    file_regression.check(contents=test_output.getvalue())
    # Test printed picture
    with open(test_image, "rb") as image:
        image_regression.check(image_data=image.read(), diff_threshold=1.0)

    # Remove image if it exists
    try:
        os.remove(test_image)
    except OSError:
        pass
