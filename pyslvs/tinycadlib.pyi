# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Union, Optional
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

def expr_parser(exprs: Sequence[Tuple[str, ...]], data_dict: Dict[str, float]) -> None:
    ...

def data_collecting(
    exprs: EStack,
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, Union[Coord, float]], int]:
    ...

def expr_solving(
    exprs: EStack,
    mapping: Dict[Union[int, Tuple[int, int]], Union[str, float]],
    vpoints: Sequence[VPoint],
    angles: Optional[Sequence[float]] = None
) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
    ...
