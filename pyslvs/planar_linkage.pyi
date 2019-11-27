# -*- coding: utf-8 -*-

from typing import Dict, Any
from numpy import ndarray
from .metaheuristics.utility import Objective

class Planar(Objective[str]):

    """This class is used to verified kinematics of the linkage mechanism."""

    def __init__(self, mech: Dict[str, Any]) -> None:
        """mech = {
            'expression': List[VPoint],
            'input': OrderedDict([((b0, d0), [start, end]), ...]),
            'placement': {pt: (x, y, r)},
            'target': {pt: [(x0, y0), (x1, y1), ...]},
            'same': {pt: match_to_pt},
            # Bound has no position data.
            'upper': List[float],
            'lower': List[float],
        }
        """
        ...

    def get_upper(self) -> ndarray:
        ...

    def get_lower(self) -> ndarray:
        ...

    def is_two_kernel(self) -> bool:
        ...

    def result(self, v: ndarray) -> str:
        """Return the last answer."""
        ...
