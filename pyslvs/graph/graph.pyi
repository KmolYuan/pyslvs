# -*- coding: utf-8 -*-

from typing import Tuple, List, Iterable, Dict
from numpy import ndarray

def link_assortment(g: Graph) -> List[int]:
    ...

def contracted_link_assortment(g: Graph) -> List[int]:
    ...

def labeled_enumerate(g: Graph) -> List[Tuple[int, Graph]]:
    ...

class Graph:

    edges: Tuple[Tuple[int, int], ...]
    vertices: Tuple[int, ...]

    def __init__(self, edges: Iterable[Tuple[int, int]]):
        """Input edges of the graph. The vertices symbols are
        positive continuously integer.
        """
        ...

    def add_edge(self, n1: int, n2: int) -> None:
        ...

    def add_vertices(self, vertices: Iterable[int]) -> None:
        ...

    def dof(self) -> int:
        ...

    def neighbors(self, n: int) -> Tuple[int, ...]:
        ...

    def degrees(self) -> Dict[int, int]:
        ...

    def degree_code(self) -> int:
        ...

    def adjacency_matrix(self) -> ndarray:
        ...

    def is_connected(self, without: int = -1) -> bool:
        ...

    def has_cut_link(self) -> bool:
        ...

    def is_degenerate(self) -> bool:
        ...

    def has_triangle(self) -> bool:
        ...

    def is_isomorphic(self, graph: Graph) -> bool:
        ...

    def is_isomorphic_vf2(self, graph: Graph) -> bool:
        ...

    def is_isomorphic_degree_code(self, graph: Graph) -> bool:
        ...

    def duplicate(self, vertices: Iterable[int], times: int) -> Graph:
        ...

    def copy(self) -> Graph:
        ...
