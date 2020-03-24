# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Optional
from .expression import VPoint

class EStack:

    def as_list(self) -> List[Tuple[str, ...]]:
        ...

def t_config(
    vpoints_: Sequence[VPoint],
    inputs: Sequence[Tuple[int, int]],
    status: Optional[Dict[int, bool]] = None
) -> EStack:
    ...
