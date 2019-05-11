# -*- coding: utf-8 -*-

from typing import (
    List,
    Tuple,
    Sequence,
    Callable,
    Optional,
)


def number_synthesis(
    nl: int,
    nj: int,
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    """Number _synthesis try-error function."""
    ...


def contracted_link(
    link_num_list: Sequence[int],
    stop_func: Optional[Callable[[], bool]] = None
) -> List[Tuple[int, ...]]:
    """Generate the contracted link assortment."""
    ...
