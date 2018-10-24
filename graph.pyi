# -*- coding: utf-8 -*-

from typing import Tuple, Sequence


class Graph:

    """NetworkX-like graph class."""

    edges: Tuple[Tuple[int, int], ...]

    def __init__(self, edges: Sequence[Tuple[int, int]]):
        ...

    def has_triangles(self) -> bool:
        ...

    def is_connected(self, with_out: int = ...) -> bool:
        ...

    def has_cut_link(self) -> bool:
        ...

    def is_isomorphic(self, graph: Graph) -> bool:
        ...

    def copy(self) -> Graph:
        ...
