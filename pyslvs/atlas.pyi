# -*- coding: utf-8 -*-

from typing import List, Sequence, Callable, Optional
from .graph import Graph

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
