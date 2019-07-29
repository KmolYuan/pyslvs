# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Dict,
    Optional,
    Union,
)
from .expression import VPoint, Coordinate

_Coord = Tuple[float, float]


class SolverSystem:

    """Sketch Solve solver."""

    def __init__(
        self,
        vpoints: Sequence[VPoint],
        inputs: Optional[Dict[Tuple[int, int], float]] = None,
        data_dict: Optional[Dict[Union[int, Tuple[int, int]], Union[Coordinate, float]]] = None
    ):
        ...

    def set_inputs(self, inputs: Dict[Tuple[int, int], float]):
        """Set input pairs."""
        ...

    def solve(self) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
        """Solve the expression."""
        ...
