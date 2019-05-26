# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    Dict,
    Union,
    Any,
)
from numpy import ndarray
from .verify import Verification


class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    def __init__(self, mech_params: Dict[str, Any]):
        """mech_params = {
            'Driver': {'pt': (x, y, r)},
            'Follower': {'pt': (x, y, r)},
            'Target': {'pt': [(x0, y0), (x1, y1), ...]},
            'constraints': [('pt', 'pt', 'pt', 'pt')],
            'Expression': str,
            'upper': ndarray[np_float32],
            'lower': ndarray[np_float32],
        }
        """
        ...

    def is_two_kernel(self) -> bool:
        ...

    def result(self, v: ndarray) -> str:
        """Return the last answer."""
        ...
