# -*- coding: utf-8 -*-

from typing import List, Tuple, Sequence, Callable, Optional

def link_synthesis(
    nl: int,
    nj: int,
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    """Return link assortment by number of links `nl` and number of joints `nj`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    ...

def contracted_link_synthesis(
    link_num_list: Sequence[int],
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    """Return contracted link assortment by link assortment `link_num_list`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    ...
