# -*- coding: utf-8 -*-

from typing import List, Tuple, Sequence, Callable, Optional

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
