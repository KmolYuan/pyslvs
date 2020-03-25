# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Callable, Optional
from .graph import Graph


def link_synthesis(
    nl: int,
    nj: int,
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    ...


def contracted_link_synthesis(
    link_num_list: Sequence[int],
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    ...

def conventional_graph(
    cg_list: List[Graph],
    c_j_list: Sequence[int],
    no_degenerate: int = 1,
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Graph]:
    ...

def contracted_graph(
    link_num: Sequence[int],
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Graph]:
    ...
