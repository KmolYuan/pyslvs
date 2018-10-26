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
    no_degenerate: bool = True,
    job_func: Optional[Callable[[List[int], int], None]] = None,
    step_func: Optional[Callable[[], None]] = None,
    stop_func: Optional[Callable[[], bool]] = None
) -> Tuple[List[Graph], float]:
    """Linkage mechanism topological function."""
    ...
