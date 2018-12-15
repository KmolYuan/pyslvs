# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    Sequence,
    Dict,
    Union,
)
from numpy import ndarray
from .verify import Verification


class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    def __init__(self, mech_params: Dict[str, Union[
        str,
        Dict[str, Tuple[float, float, float]],
        Dict[str, Sequence[Tuple[float, float]]],
        Sequence[Tuple[str, str, str, str]],
        ndarray,
    ]]):
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

    def result(self, v: ndarray) -> Dict[str, Union[Tuple[float, float], float]]:
        """Return the last answer."""
        ...
