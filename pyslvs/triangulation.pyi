# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Optional
from .expression import VPoint

class ExpressionStack:

    """The stack of Python wrapper."""

    def as_list(self) -> List[Tuple[str, ...]]:
        """Turn into normal list."""
        ...

    def __repr__(self) -> str:
        ...

def vpoints_configure(
    vpoints_: Sequence[VPoint],
    inputs: Sequence[Tuple[int, int]],
    status: Optional[Dict[int, bool]] = None
) -> ExpressionStack:
    """Auto configuration algorithm."""
    ...
