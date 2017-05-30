r"""
Parsers read output files from MD simulation packages and create 
SimulationData objects with their contents.
"""


class Parser(object):
    r"""
    Parser base class
    """
    def get_simulation_data(self):
        raise NotImplementedError
