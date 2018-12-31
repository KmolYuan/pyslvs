# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Dict,
    Union,
    Optional,
)
from .pmks import VPoint

TuplePoint = Tuple[float, float]


class Coordinate:

    """A class to store the coordinate."""

    x: float
    y: float

    def __init__(self, x: float, y: float):
        ...

    def distance(self, p: Coordinate) -> float:
        """Distance."""
        ...

    def is_nan(self) -> bool:
        """Test this coordinate is a error-occurred answer."""
        ...

    def __repr__(self) -> str:
        ...

def plap(
    c1: Coordinate,
    d0: float,
    a0: float,
    c2: Optional[Coordinate] = None,
    inverse: bool = False
) -> TuplePoint:
    """Point on circle by angle."""
    ...

def pllp(
    c1: Coordinate,
    d0: float,
    d1: float,
    c2: Coordinate,
    inverse: bool = False
) -> TuplePoint:
    """Two intersection points of two circles."""
    ...

def plpp(
    c1: Coordinate,
    d0: float,
    c2: Coordinate,
    c3: Coordinate,
    inverse: bool = False
) -> TuplePoint:
    """Two intersection points of a line and a circle."""
    ...

def pxy(c1: Coordinate, x: float, y: float) -> TuplePoint:
    """Using relative cartesian coordinate to get solution."""
    ...

def vpoint_dof(vpoints: Sequence[VPoint]) -> int:
    """Degree of freedoms calculate from PMKS expressions."""
    ...

def expr_parser(exprs: Sequence[Tuple[str, ...]], data_dict: Dict[str, float]):
    """Update data."""
    ...

def expr_solving(
    exprs: Sequence[Tuple[str, ...]],
    mapping: Dict[int, str],
    vpoints: Sequence[VPoint],
    angles: Sequence[float] = None
) -> List[Union[TuplePoint, Tuple[TuplePoint, TuplePoint]]]:
    """Solving function."""
    ...

def data_collecting(
    exprs: Sequence[Tuple[str, ...]],
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, float], int]:
    """Data collecting process."""
    ...
