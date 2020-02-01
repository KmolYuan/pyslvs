# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Optional
from .expression import VPoint

class ExpressionStack:

    """Triangle solution stack, generated from [`vpoints_configure`](#vpoints_configure).
    It is pointless to call the constructor.
    """

    def as_list(self) -> List[Tuple[str, ...]]:
        """Copy the dataset as list object."""
        ...

    def __repr__(self) -> str:
        ...

def vpoints_configure(
    vpoints_: Sequence[VPoint],
    inputs: Sequence[Tuple[int, int]],
    status: Optional[Dict[int, bool]] = None
) -> ExpressionStack:
    """Generate the Triangle solution stack by mechanism expression `vpoints_`.

    The argument `inputs` is a list of input pairs.
    The argument `status` will track the configuration of each point, which is optional.
    """
    ...
