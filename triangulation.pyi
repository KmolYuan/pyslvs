# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Sequence,
    Dict,
    Optional,
)
from .pmks import VPoint


def vpoints_configure(
    vpoints_: Sequence[VPoint],
    inputs: Sequence[Tuple[int, int]],
    status: Optional[Dict[int, bool]] = None
) -> List[Tuple[str, ...]]:
    """Auto configuration algorithm."""
    ...
