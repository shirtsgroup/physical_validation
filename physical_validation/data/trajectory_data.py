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
Data structures carrying simulation data.
"""
import warnings

import numpy as np

from ..util import error as pv_error
from ..util.util import array_equal_shape_and_close


class Box(object):
    def __init__(self, box=None):
        pass

    def get(self, key):
        raise NotImplementedError

    def __getitem__(self, key):
        raise NotImplementedError

    def set(self, key, value):
        raise NotImplementedError

    def __setitem__(self, key, value):
        raise NotImplementedError

    @property
    def volume(self):
        raise NotImplementedError

    @property
    def box(self):
        raise NotImplementedError

    @box.setter
    def box(self, b):
        raise NotImplementedError

    def gather(self, positions, bonds, molec_idx):
        raise NotImplementedError


class RectangularBox(Box):
    def __init__(self, box=None):
        Box.__init__(self)
        self.__box = None
        self.__nframes = 0

        if box is not None:
            self.box = box

        self.__getters = {
            "box": RectangularBox.box.__get__,
            "volume": RectangularBox.volume.__get__,
        }
        self.__setters = {"box": RectangularBox.box.__set__}

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.__getters:
            raise KeyError
        return self.__getters[key](self)

    def set(self, key, value):
        self[key] = value

    def __setitem__(self, key, value):
        if key not in self.__setters:
            raise KeyError
        self.__setters[key](self, value)

    @property
    def volume(self):
        return np.prod(self.box, axis=1)

    @property
    def box(self):
        return self.__box

    @box.setter
    def box(self, b):
        b = np.array(b)
        assert 0 < b.ndim < 3
        if b.ndim == 1:
            assert b.size == 3
            self.__box = b
            self.__nframes = 1
        elif b.ndim == 2:
            assert b.shape[1] == 3
            self.__box = b
            self.__nframes = b.shape[0]

    def gather(self, positions, bonds, molec_idx):
        bonds = np.array(bonds)
        if bonds.size == 0:
            return positions
        positions = np.array(positions)
        assert 1 < positions.ndim < 4
        if positions.ndim == 2:
            nframes = 1
            positions = np.array([positions])
        else:
            nframes = positions.shape[0]
        if self.__nframes != 1:
            assert self.__nframes == nframes
        for f in range(nframes):
            p = positions[f]
            if self.__nframes > 1:
                box = self.box[f]
            else:
                box = self.box[0]
            assert len(bonds) == len(molec_idx)
            for mbonds, idx in zip(bonds, molec_idx):
                for b in mbonds:
                    a1 = idx + b[0]
                    a2 = idx + b[1]
                    p[a2] += np.round((p[a1] - p[a2]) / box) * box
            positions[f] = p
        return positions


class TrajectoryData(object):
    r"""TrajectoryData: The position and velocity trajectory along the simulation

    The full trajectory is needed to calculate the equipartition of the kinetic energy.
    As they are used in connection, the position and velocity trajectories are expected
    to have the same shape and number of frames.

    The position and velocity trajectories can be accessed either using the getters
    of an object, as in

        * trajectory.position
        * trajectory.velocity

    or using the key notation, as in

        * trajectory['position']
        * trajectory['velocity']

    """

    @staticmethod
    def trajectories():
        return ("position", "velocity")

    def __init__(self, position=None, velocity=None):
        self.__position = None
        self.__velocity = None
        self.__nframes = 0

        if position is not None:
            self.position = position
        if velocity is not None:
            self.velocity = velocity

        self.__getters = {
            "position": TrajectoryData.position.__get__,
            "velocity": TrajectoryData.velocity.__get__,
        }

        self.__setters = {
            "position": TrajectoryData.position.__set__,
            "velocity": TrajectoryData.velocity.__set__,
        }

    def get(self, key):
        return self[key]

    def __getitem__(self, key):
        if key not in self.trajectories():
            raise KeyError
        return self.__getters[key](self)

    def set(self, key, value):
        self[key] = value

    def __setitem__(self, key, value):
        if key not in self.trajectories():
            raise KeyError
        self.__setters[key](self, value)

    @property
    def position(self):
        """Get position"""
        return self.__position

    @position.setter
    def position(self, pos):
        """Set position"""
        pos = np.array(pos)
        if pos.ndim == 2:
            # create 3-dimensional array
            pos = np.array([pos])
        if pos.ndim != 3:
            warnings.warn("Expected 2- or 3-dimensional array.")
        if self.__nframes == 0 and self.__velocity is None:
            self.__nframes = pos.shape[0]
        elif self.__nframes != pos.shape[0]:
            raise pv_error.InputError(
                ["pos"], "Expected equal number of frames as in velocity trajectory."
            )
        self.__position = pos

    @property
    def velocity(self):
        """Get velocity"""
        return self.__velocity

    @velocity.setter
    def velocity(self, vel):
        """Set velocity"""
        vel = np.array(vel)
        if vel.ndim == 2:
            # create 3-dimensional array
            vel = np.array([vel])
        if vel.ndim != 3:
            warnings.warn("Expected 2- or 3-dimensional array.")
        if self.__nframes == 0 and self.__position is None:
            self.__nframes = vel.shape[0]
        elif self.__nframes != vel.shape[0]:
            raise pv_error.InputError(
                ["vel"], "Expected equal number of frames as in position trajectory."
            )
        self.__velocity = vel

    @property
    def nframes(self):
        """Get number of frames"""
        return self.__nframes

    def __eq__(self, other):
        if type(other) is not type(self):
            return False

        return (
            array_equal_shape_and_close(self.__position, other.__position)
            and array_equal_shape_and_close(self.__velocity, other.__velocity)
            and self.__nframes == other.__nframes
        )
