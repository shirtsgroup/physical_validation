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
Unit and regression test for the physical_validation package.
"""


def test_physical_validation_imported():
    """Sample test, will always pass so long as import statement worked"""

    import physical_validation

    assert physical_validation is not None
    assert physical_validation.data is not None
    assert physical_validation.ensemble is not None
    assert physical_validation.integrator is not None
    assert physical_validation.kinetic_energy is not None
    assert physical_validation.util is not None
