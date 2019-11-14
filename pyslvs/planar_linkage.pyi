# -*- coding: utf-8 -*-

from typing import Dict, Any
from numpy import ndarray
from .Adesign.verify import Verification

class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    def __init__(self, mech: Dict[str, Any]) -> None:
        """mech = {
            'expression': List[VPoint],
            'input': {(b0, d0): (start, end), ...},
            'placement': {pt: (x, y, r)},
            'target': {pt: [(x0, y0), (x1, y1), ...]},
            'same': {pt: match_to_pt},
            # Bound has no position data.
            'upper': List[float],
            'lower': List[float],
        }
        """
        ...

    def is_two_kernel(self) -> bool:
        ...

    def result(self, v: ndarray) -> str:
        """Return the last answer."""
        ...
