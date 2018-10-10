# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Dict,
    Mapping,
    Optional,
)
from pmks import VPoint


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

def PLAP(
    A: Coordinate,
    L0: float,
    a0: float,
    B: Optional[Coordinate],
    inverse: bool = False
) -> Tuple[float, float]:
    """Point on circle by angle."""
    ...

def PLLP(
    A: Coordinate,
    L0: float,
    L1: float,
    B: Coordinate,
    inverse: bool = False
) -> Tuple[float, float]:
    """Two intersection points of two circles."""
    ...

def PLPP(
    A: Coordinate,
    L0: float,
    B: Coordinate,
    C: Coordinate,
    inverse: bool = False
) -> Tuple[float, float]:
    """Two intersection points of a line and a circle."""
    ...

def PXY(A: Coordinate, x: float, y: float) -> Tuple[float, float]:
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
    angles: Mapping[int, float] = None
) -> List[Tuple[float, float]]:
    """Solving function."""
    ...

def data_collecting(
    exprs: Sequence[Tuple[str, ...]],
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, float], int]:
    """Data collecting process."""
    ...
