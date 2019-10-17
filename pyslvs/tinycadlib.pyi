# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Union, Optional
from .triangulation import ExpressionStack
from .expression import VPoint, Coordinate

TuplePoint = Tuple[float, float]

def plap(
    c1: Coordinate,
    d0: float,
    a0: float,
    c2: Optional[Coordinate] = None,
    inverse: bool = False
) -> Coordinate:
    """Point on circle by angle."""
    ...

def pllp(
    c1: Coordinate,
    d0: float,
    d1: float,
    c2: Coordinate,
    inverse: bool = False
) -> Coordinate:
    """Two intersection points of two circles."""
    ...

def plpp(
    c1: Coordinate,
    d0: float,
    c2: Coordinate,
    c3: Coordinate,
    inverse: bool = False
) -> Coordinate:
    """Two intersection points of a line and a circle."""
    ...

def pxy(c1: Coordinate, x: float, y: float) -> Coordinate:
    """Using relative cartesian coordinate to get solution."""
    ...

def vpoint_dof(vpoints: Sequence[VPoint]) -> int:
    """Degree of freedoms calculate from PMKS expressions."""
    ...

def expr_parser(exprs: Sequence[Tuple[str, ...]], data_dict: Dict[str, float]) -> None:
    """Update data."""
    ...

def data_collecting(
    exprs: ExpressionStack,
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, Union[Coordinate, float]], int]:
    """Data collecting process."""
    ...

def expr_solving(
    exprs: ExpressionStack,
    mapping: Dict[Union[int, Tuple[int, int]], Union[str, float]],
    vpoints: Sequence[VPoint],
    angles: Optional[Sequence[float]] = None
) -> List[Union[TuplePoint, Tuple[TuplePoint, TuplePoint]]]:
    """Solving function."""
    ...
