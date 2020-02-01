# -*- coding: utf-8 -*-

from typing import List, Sequence, Callable, Optional
from .graph import Graph

def conventional_graph(
    cg_list: List[Graph],
    c_j_list: Sequence[int],
    no_degenerate: int = 1,
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Graph]:
    """Generate conventional graphs by contracted graphs `cg_list` and
    contracted link assortment `c_j_list`.

    The degenerate setting `no_degenerate` has following option:

    + `0`: No degenerate.
    + `1`: Only degenerate.
    + Else: All graphs.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    ...

def contracted_graph(
    link_num: Sequence[int],
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Graph]:
    """Generate contracted graphs by link assortment `link_num`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    ...
