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
    job_func: Optional[Callable[[List[int], int], None]] = None,
    step_func: Optional[Callable[[], None]] = None,
    stop_func: Optional[Callable[[], bool]] = None
) -> Tuple[List[Graph], float]:
    """Linkage mechanism topological function."""
    ...
