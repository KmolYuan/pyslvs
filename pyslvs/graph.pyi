# -*- coding: utf-8 -*-

from typing import Tuple, List, Iterable, Dict
from numpy import ndarray

def link_assortment(g: Graph) -> List[int]:
    """Return link assortment of the graph."""
    ...

def contracted_link_assortment(g: Graph) -> List[int]:
    """Return contracted link assortment of the graph."""
    ...

def labeled_enumerate(g: Graph) -> List[Tuple[int, Graph]]:
    """Enumerate each node with labeled except isomorphism."""
    ...

class Graph:

    """NetworkX-like graph class."""

    edges: Tuple[Tuple[int, int], ...]
    vertices: Tuple[int, ...]

    def __init__(self, edges: Iterable[Tuple[int, int]]) -> None: ...

    def add_edge(self, n1: int, n2: int) -> None:
        """Add two vertices for an edge."""
        ...

    def add_vertices(self, vertices: Iterable[int]) -> None:
        """Add vertices from a tuple."""
        ...

    def dof(self) -> int:
        """Return degrees of freedom."""
        ...

    def neighbors(self, n: int) -> Tuple[int, ...]:
        """Neighbors except the node."""
        ...

    def degrees(self) -> Dict[int, int]:
        """Neighbors except the node."""
        ...

    def degree_code(self) -> int:
        """Return degree code of the graph."""
        ...

    def adjacency_matrix(self) -> ndarray:
        """Represent as adjacency matrix."""
        ...

    def is_connected(self, without: int = -1) -> bool:
        """Return True if the graph is not isolated."""
        ...

    def has_cut_link(self) -> bool:
        """Return True if the graph has any cut links."""
        ...

    def is_degenerate(self) -> bool:
        """Return True if this kinematic chain is degenerate.

        + Prue all multiple contracted links recursively.
        + Check the DOF of sub-graph if it is lower then zero.
        """
        ...

    def has_triangle(self) -> bool:
        """Return True if the graph has triangle."""
        ...

    def is_isomorphic(self, graph: Graph) -> bool:
        """Return True if two graphs is isomorphic."""
        ...

    def is_isomorphic_vf2(self, graph: Graph) -> bool:
        """Compare isomorphism by VF2 algorithm."""
        ...

    def is_isomorphic_degree_code(self, graph: Graph) -> bool:
        """Compare isomorphism by degree code algorithm."""
        ...

    def duplicate(self, vertices: Iterable[int], times: int) -> Graph:
        """Make graph duplicate by specific vertices. Return a new graph."""
        ...

    def copy(self) -> Graph:
        """Copy the graph."""
        ...
