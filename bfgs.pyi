# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Dict,
    Optional,
    Union,
)
from .pmks import VPoint

def vpoint_solving(
    vpoints: Sequence[VPoint],
    inputs: Optional[Dict[Tuple[int, int], float]] = None,
    data_dict: Optional[Dict[Union[int, Tuple[int, int]], float]] = None
) -> List[Tuple[float, float]]:
    """Solving function from vpoint list."""
    ...
