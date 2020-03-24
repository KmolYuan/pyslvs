# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Union, Optional
from .triangulation import EStack
from .expression import VPoint, Coordinate

TuplePoint = Tuple[float, float]

def plap(
    c1: Coordinate,
    d0: float,
    a0: float,
    c2: Optional[Coordinate] = None,
    inverse: bool = False
) -> Coordinate:
    ...

def pllp(
    c1: Coordinate,
    d0: float,
    d1: float,
    c2: Coordinate,
    inverse: bool = False
) -> Coordinate:
    ...

def plpp(
    c1: Coordinate,
    d0: float,
    c2: Coordinate,
    c3: Coordinate,
    inverse: bool = False
) -> Coordinate:
    ...

def pxy(c1: Coordinate, x: float, y: float) -> Coordinate:
    ...

def vpoint_dof(vpoints: Sequence[VPoint]) -> int:
    ...

def expr_parser(exprs: Sequence[Tuple[str, ...]], data_dict: Dict[str, float]) -> None:
    ...

def data_collecting(
    exprs: EStack,
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, Union[Coordinate, float]], int]:
    ...

def expr_solving(
    exprs: EStack,
    mapping: Dict[Union[int, Tuple[int, int]], Union[str, float]],
    vpoints: Sequence[VPoint],
    angles: Optional[Sequence[float]] = None
) -> List[Union[TuplePoint, Tuple[TuplePoint, TuplePoint]]]:
    ...
