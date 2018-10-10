# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    Sequence,
    Callable,
    Optional,
)
from graph import Graph


def topo(
    link_num: Sequence[int],
    degenerate: bool = True,
    job_func: Optional[Callable[[str, int], None]] = None,
    stop_func: Optional[Callable[[], None]] = None
) -> Tuple[Graph]:
    """Linkage mechanism topological function."""
    ...
