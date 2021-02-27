###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (c) 2017-2021 University of Colorado Boulder               #
#              (c) 2012      The University of Virginia                   #
#                                                                         #
###########################################################################
r"""
Module containing the custom exception classes for the physical_validation
package.
"""


class PhysicalValidationError(Exception):
    r"""Base class for exceptions in the physical_validation module."""
    pass


class InputError(PhysicalValidationError):
    r"""Exception raised for input errors"""

    def __init__(self, argument, message):
        r"""

        Parameters
        ----------
        argument : string or list of strings
        message : string
        """
        self.argument = argument
        self.message = message


class ParserValueNotSetError(PhysicalValidationError):
    r"""
    Exception raised if a requested data value
    was not set by the user previously
    """

    def __init__(self, message):
        r"""

        Parameters
        ----------
        message : string
        """
        self.message = message


class FileFormatError(PhysicalValidationError):
    r"""Exception raised for files not following expected format"""

    def __init__(self, argument, message):
        r"""

        Parameters
        ----------
        argument : string or list of strings
        message : string
        """
        self.argument = argument
        self.message = message
