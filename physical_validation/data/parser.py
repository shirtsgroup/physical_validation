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
Parsers read output files from MD simulation packages and create
SimulationData objects with their contents.
"""
from . import SimulationData


class Parser(object):
    r"""
    Parser base class
    """

    def get_simulation_data(self) -> SimulationData:
        raise NotImplementedError
