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
r"""
Regression tests running the ensemble checks and making sure that the
output is unchanged between code versions.

Output changes might be intended. If you are sure a change is wanted,
you can update the reference data using
    $ pytest test_ensemble_regression.py --force-regen
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
from typing import List, Optional, Tuple

import pytest

import physical_validation as pv

from .test_system_database import database


def run_ensemble_check(
    simulation_data_1: pv.data.SimulationData,
    simulation_data_2: pv.data.SimulationData,
    image_filename: Optional[str] = None,
) -> List[float]:
    r"""
    Runs the ensemble check on the provided simulation data objects.

    Parameters
    ----------
    simulation_data_1
        First simulation data object
    simulation_data_2
        Second simulation data object
    image_filename
        Plot distributions to `filename`.
        Default: None, no plotting to file.


    Returns
    -------
    Forwards the output from the ensemble check.
    """
    print("## Validating ensemble")
    # Run ensemble check
    quantiles = pv.ensemble.check(
        simulation_data_1,
        simulation_data_2,
        verbosity=2,
        filename=image_filename,
        bs_seed=1,
    )

    # Print result
    if len(quantiles) == 1:
        q_str = "{:.1f}".format(quantiles[0])
    else:
        q_str = "("
        for q in quantiles:
            q_str += "{:.1f}, ".format(q)
        q_str = q_str[:-2] + ")"
    print("Calculated slope is " + q_str + " quantiles from the true slope")

    return quantiles


def ensemble_nvt_flat_file(image_filename: Optional[str] = None) -> List[float]:
    r"""
    Read NVT flat file data and launch the ensemble checks.

    Parameters
    ----------
    image_filename
        Plot distributions to `filename`.

    Returns
    -------
    Forwards the output from the ensemble check.
    """
    system_name = "Water900"
    print("### Regression test of NVT ensemble using flat files")
    print("### System: " + system_name)

    parser = pv.data.FlatfileParser()
    system = database.system(system_name)
    print("## Reading low temperature result")
    simulation_data_low = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble("NVT-low"),
        system=system.system_data,
        kinetic_ene_file=system.observable_flat_file("NVT-low", "kinetic_energy"),
        potential_ene_file=system.observable_flat_file("NVT-low", "potential_energy"),
        total_ene_file=system.observable_flat_file("NVT-low", "total_energy"),
        volume_file=system.observable_flat_file("NVT-low", "volume"),
    )
    print("## Reading high temperature result")
    simulation_data_high = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble("NVT-high"),
        system=system.system_data,
        kinetic_ene_file=system.observable_flat_file("NVT-high", "kinetic_energy"),
        potential_ene_file=system.observable_flat_file("NVT-high", "potential_energy"),
        total_ene_file=system.observable_flat_file("NVT-high", "total_energy"),
        volume_file=system.observable_flat_file("NVT-high", "volume"),
    )

    return run_ensemble_check(
        simulation_data_1=simulation_data_low,
        simulation_data_2=simulation_data_high,
        image_filename=image_filename,
    )


def ensemble_nvt_numpy_arrays(image_filename: Optional[str] = None) -> List[float]:
    r"""
    Create NVT data in numpy arrays and launch the ensemble checks.

    Parameters
    ----------
    image_filename
        Plot distributions to `filename`.

    Returns
    -------
    Forwards the output from the ensemble check.
    """
    system_name = "Water900"
    print("### Regression test of NVT ensemble using numpy arrays")
    print("### System: " + system_name)

    system = database.system(system_name)
    print("## Creating low temperature result")
    observables = pv.data.ObservableData(
        kinetic_energy=system.observable_as_array("NVT-low", "kinetic_energy"),
        potential_energy=system.observable_as_array("NVT-low", "potential_energy"),
        total_energy=system.observable_as_array("NVT-low", "total_energy"),
        volume=system.observable_as_array("NVT-low", "volume"),
    )
    simulation_data_low = pv.data.SimulationData(
        units=system.units,
        ensemble=system.ensemble("NVT-low"),
        system=system.system_data,
        observables=observables,
    )
    print("## Creating high temperature result")
    observables = pv.data.ObservableData(
        kinetic_energy=system.observable_as_array("NVT-high", "kinetic_energy"),
        potential_energy=system.observable_as_array("NVT-high", "potential_energy"),
        total_energy=system.observable_as_array("NVT-high", "total_energy"),
        volume=system.observable_as_array("NVT-high", "volume"),
    )
    simulation_data_high = pv.data.SimulationData(
        units=system.units,
        ensemble=system.ensemble("NVT-high"),
        system=system.system_data,
        observables=observables,
    )

    return run_ensemble_check(
        simulation_data_1=simulation_data_low,
        simulation_data_2=simulation_data_high,
        image_filename=image_filename,
    )


def get_npt_simulation_ids(identifier: str) -> Tuple[str, str]:
    r"""
    Return the appropriate simulation ids given
    Parameters
    ----------
    identifier

    Returns
    -------

    """
    if identifier == "Temperature only":
        return "NPT-lowT-lowP", "NPT-highT-lowP"
    elif identifier == "Pressure only":
        return "NPT-lowT-lowP", "NPT-lowT-highP"
    elif identifier == "Temperature and pressure":
        return "NPT-lowT-lowP", "NPT-highT-highP"
    else:
        raise KeyError("identifier")


def ensemble_npt_flat_file(
    test_type: str, image_filename: Optional[str] = None
) -> List[float]:
    r"""
    Read NPT flat file data and launch the ensemble checks.

    Parameters
    ----------
    test_type
        The identifier of the type of simulation results we're reading.
    image_filename
        Plot distributions to `filename`.

    Returns
    -------
    Forwards the output from the ensemble check.
    """
    system_name = "Water900"
    print("### Regression test of NPT ensemble using flat files")
    print("### System: " + system_name)

    id_low, id_high = get_npt_simulation_ids(test_type)

    parser = pv.data.FlatfileParser()
    system = database.system(system_name)
    print("## Reading low temperature result")
    simulation_data_low = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble(id_low),
        system=system.system_data,
        kinetic_ene_file=system.observable_flat_file(id_low, "kinetic_energy"),
        potential_ene_file=system.observable_flat_file(id_low, "potential_energy"),
        total_ene_file=system.observable_flat_file(id_low, "total_energy"),
        volume_file=system.observable_flat_file(id_low, "volume"),
    )
    print("## Reading high temperature result")
    simulation_data_high = parser.get_simulation_data(
        units=system.units,
        ensemble=system.ensemble(id_high),
        system=system.system_data,
        kinetic_ene_file=system.observable_flat_file(id_high, "kinetic_energy"),
        potential_ene_file=system.observable_flat_file(id_high, "potential_energy"),
        total_ene_file=system.observable_flat_file(id_high, "total_energy"),
        volume_file=system.observable_flat_file(id_high, "volume"),
    )

    return run_ensemble_check(
        simulation_data_1=simulation_data_low,
        simulation_data_2=simulation_data_high,
        image_filename=image_filename,
    )


def ensemble_npt_numpy_arrays(
    test_type: str, image_filename: Optional[str] = None
) -> List[float]:
    r"""
    Create NPT data in numpy arrays and launch the ensemble checks.

    Parameters
    ----------
    test_type
        The identifier of the type of simulation results we're reading.
    image_filename
        Plot distributions to `filename`.

    Returns
    -------
    Forwards the output from the ensemble check.
    """
    system_name = "Water900"
    print("### Regression test of NPT ensemble using numpy arrays")
    print("### System: " + system_name)

    id_low, id_high = get_npt_simulation_ids(test_type)

    system = database.system(system_name)
    print("## Creating low temperature result")
    observables = pv.data.ObservableData(
        kinetic_energy=system.observable_as_array(id_low, "kinetic_energy"),
        potential_energy=system.observable_as_array(id_low, "potential_energy"),
        total_energy=system.observable_as_array(id_low, "total_energy"),
        volume=system.observable_as_array(id_low, "volume"),
    )
    simulation_data_low = pv.data.SimulationData(
        units=system.units,
        ensemble=system.ensemble(id_low),
        system=system.system_data,
        observables=observables,
    )
    print("## Creating high temperature result")
    observables = pv.data.ObservableData(
        kinetic_energy=system.observable_as_array(id_high, "kinetic_energy"),
        potential_energy=system.observable_as_array(id_high, "potential_energy"),
        total_energy=system.observable_as_array(id_high, "total_energy"),
        volume=system.observable_as_array(id_high, "volume"),
    )
    simulation_data_high = pv.data.SimulationData(
        units=system.units,
        ensemble=system.ensemble(id_high),
        system=system.system_data,
        observables=observables,
    )

    return run_ensemble_check(
        simulation_data_1=simulation_data_low,
        simulation_data_2=simulation_data_high,
        image_filename=image_filename,
    )


@pytest.mark.parametrize("input_source", ["flat file", "numpy array"])
def test_ensemble_regression_nvt(
    data_regression, file_regression, image_regression, input_source: str
) -> None:
    r"""
    Regression test running NVT ensemble checks.

    Parameters
    ----------
    data_regression
        Regression test fixture testing python dicts
    file_regression
        Regression test fixture testing text files
    image_regression
        Regression test fixture testing images
    input_source
        Whether we're using the flat file parsers or numpy arrays
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
        if input_source == "flat file":
            result = ensemble_nvt_flat_file(image_filename=test_image)
        elif input_source == "numpy array":
            result = ensemble_nvt_numpy_arrays(image_filename=test_image)
        else:
            raise NotImplementedError("Unknown input source " + input_source)

    # Test returned value (regression is only checking dicts of strings)
    result_dict = {
        n: "{:.6f}".format(list_entry) for n, list_entry in enumerate(result)
    }
    data_regression.check(result_dict)
    # Test printed output
    file_regression.check(contents=test_output.getvalue())
    # Test printed picture
    with open(test_image, "rb") as image:
        image_regression.check(image_data=image.read(), diff_threshold=0.5)

    # Clean up
    try:
        os.remove(test_image)
    except OSError:
        pass


@pytest.mark.parametrize("input_source", ["flat file", "numpy array"])
@pytest.mark.parametrize(
    "test_type", ["Temperature only", "Pressure only", "Temperature and pressure"]
)
def test_ensemble_regression_npt(
    data_regression,
    file_regression,
    image_regression,
    input_source: str,
    test_type: str,
) -> None:
    r"""
    Regression test running NVT ensemble checks.

    Parameters
    ----------
    data_regression
        Regression test fixture testing python dicts
    file_regression
        Regression test fixture testing text files
    image_regression
        Regression test fixture testing images
    input_source
        Whether we're using the flat file parsers or numpy arrays
    test_type
        Whether we're testing results at different temperatures, different
        pressures, or both different temperatures and pressures
    """
    test_output = StringIO()
    is_2d = test_type == "Temperature and pressure"
    # no plotting in 2D
    test_image = "test_plot.png" if not is_2d else None

    # Remove image if it exists
    if test_image is not None:
        try:
            os.remove(test_image)
        except OSError:
            pass

    # Redirect stdout into string which we can test
    with redirect_stdout(test_output):
        if input_source == "flat file":
            result = ensemble_npt_flat_file(
                test_type=test_type, image_filename=test_image
            )
        elif input_source == "numpy array":
            result = ensemble_npt_numpy_arrays(
                test_type=test_type, image_filename=test_image
            )
        else:
            raise NotImplementedError("Unknown input source " + input_source)

    # Test returned value (regression is only checking dicts of strings)
    result_dict = {
        n: "{:.6f}".format(list_entry) for n, list_entry in enumerate(result)
    }
    data_regression.check(result_dict)
    # Test printed output
    file_regression.check(contents=test_output.getvalue())
    # Test printed picture
    if test_image is not None:
        with open(test_image, "rb") as image:
            image_regression.check(image_data=image.read(), diff_threshold=0.5)

    # Clean up
    if test_image is not None:
        try:
            os.remove(test_image)
        except OSError:
            pass
