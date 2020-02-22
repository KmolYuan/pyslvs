# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Optional
from .expression import VPoint

class ExpressionStack:

    def as_list(self) -> List[Tuple[str, ...]]:
        ...

def vpoints_configure(
    vpoints_: Sequence[VPoint],
    inputs: Sequence[Tuple[int, int]],
    status: Optional[Dict[int, bool]] = None
) -> ExpressionStack:
    ...
