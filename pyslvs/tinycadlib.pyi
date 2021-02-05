# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Mapping, Union, Optional
from numpy import ndarray
from .topo_config import EStack
from .expression import VPoint, Coord

_Coord = Tuple[float, float]

def pxy(c1: Coord, x: float, y: float) -> Coord:
    ...

def plap(
    c1: Coord,
    d0: float,
    a0: float,
    c2: Optional[Coord] = None,
    inverse: bool = False
) -> Coord:
    ...

def ppp(c1: Coord, c2: Coord, c3: Coord) -> Coord:
    ...

def pllp(c1: Coord, d0: float, d1: float, c2: Coord,
         inverse: bool = False) -> Coord:
    ...

def plpp(c1: Coord, d0: float, c2: Coord, c3: Coord,
         inverse: bool = False) -> Coord:
    ...

def palp(c1: Coord, a0: float, d0: float, c2: Coord, inverse: bool = False) -> Coord:
    ...

def vpoint_dof(vpoints: Sequence[VPoint]) -> int:
    ...

def expr_solving(
    exprs: EStack,
    vpoints: Sequence[VPoint],
    angles: Optional[Mapping[Tuple[int, int], float]] = None
) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
    ...

def uniform_four_bar(ml: float, n: int) -> ndarray:
    ...

def uniform_path(v: ndarray, n: int) -> ndarray:
    ...

def uniform_expr(v: ndarray) -> List[VPoint]:
    ...
