# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Callable,
    Optional,
)
from .graph import Graph


def topo(
    link_num: Sequence[int],
    c_j_list: Sequence[int],
    no_degenerate: int = 1,
    stop_func: Optional[Callable[[], bool]] = None
) -> Tuple[List[Graph], float]:
    """Linkage mechanism topological function."""
    ...
