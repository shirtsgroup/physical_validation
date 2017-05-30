

class PhysicalValidationError(Exception):
    r"""Base class for exceptions in the physicalvalidation module."""
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
