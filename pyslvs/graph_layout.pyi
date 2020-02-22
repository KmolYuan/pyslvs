# -*- coding: utf-8 -*-

from typing import Tuple, Dict
from .graph import Graph

def external_loop_layout(
    graph: Graph,
    node_mode: bool,
    scale: float = 1.
) -> Dict[int, Tuple[float, float]]:
    ...
