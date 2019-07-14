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


def vpoint_solving(
    vpoints: Sequence[VPoint],
    inputs: Optional[Dict[Tuple[int, int], float]] = None,
    data_dict: Optional[Dict[Union[int, Tuple[int, int]], Union[Coordinate, float]]] = None
) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
    """Solving function from vpoint list."""
    ...
