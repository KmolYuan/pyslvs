# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Iterable,
)


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
    nodes: Tuple[int, ...]

    def __init__(self, edges: Iterable[Tuple[int, int]]): ...

    def add_edge(self, n1: int, n2: int):
        """Add two nodes for an edge."""
        ...

    def add_nodes(self, nodes: Iterable[int]):
        """Add nodes from a tuple."""
        ...

    def dof(self) -> int:
        """Return degrees of freedom."""
        ...

    def neighbors(self, n: int) -> Tuple[int, ...]:
        """Neighbors except the node."""
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

    def is_isomorphic(self, graph: Graph) -> bool:
        """Return True if two graphs is isomorphic."""
        ...

    def duplicate(self, nodes: Iterable[int], times: int) -> Graph:
        """Make graph duplicate by specific nodes. Return a new graph."""
        ...

    def copy(self) -> Graph:
        """Copy the graph."""
        ...
