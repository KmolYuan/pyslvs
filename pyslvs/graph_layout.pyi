# -*- coding: utf-8 -*-

from typing import Tuple, Dict
from .graph import Graph

def external_loop_layout(
    graph: Graph,
    node_mode: bool,
    scale: float = 1.
) -> Dict[int, Tuple[float, float]]:
    """Layout position decided by outer loop (max cycle).

    Return the layout position decided by external loop.
    Argument `node_mode` will transform edges into vertices.
    Argument `scale` will resize the position by scale factor.
    """
    ...
