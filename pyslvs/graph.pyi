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

    """The undirected graph class, support multigraph."""

    edges: Tuple[Tuple[int, int], ...]
    vertices: Tuple[int, ...]

    def __init__(self, edges: Iterable[Tuple[int, int]]) -> None:
        """Input edges of the graph. The vertices symbols are
        positive continuously integer.
        """
        ...

    def add_edge(self, n1: int, n2: int) -> None:
        """Add edge `n1` to `n2`."""
        ...

    def add_vertices(self, vertices: Iterable[int]) -> None:
        """Add vertices from iterable object `vertices`."""
        ...

    def dof(self) -> int:
        """Return DOF of the graph.

        !!! note
            DOF is the Degree of Freedoms to a mechanism.

            In the [Graph] objects, all vertices will assumed as revolute joints (1 DOF).

            $$
            F = 3(N_L - 1) - 2N_J
            $$
        """
        ...

    def neighbors(self, n: int) -> Tuple[int, ...]:
        """Return the neighbors of the vertex `n`."""
        ...

    def degrees(self) -> Dict[int, int]:
        """Return the degrees of each vertex."""
        ...

    def degree_code(self) -> int:
        """Generate a degree code.

        With a sorted vertices mapping by the degrees of each vertex,
        regenerate a new adjacency matrix.
        A binary code can be found by concatenating the upper right elements.
        The degree code is the maximum value of the permutation.
        """
        ...

    def adjacency_matrix(self) -> ndarray:
        """Generate a adjacency matrix.

        Assume the matrix $A[i, j] = A[j, i]$.
        Where $A[i, j] = 1$ if edge `(i, j)` exist.
        """
        ...

    def is_connected(self, without: int = -1) -> bool:
        """Return `True` if the graph is connected.
        Set the argument `without` to ignore one vertex.
        """
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
        """Return True if two graphs is isomorphic.

        Default is using VF2 algorithm.
        """
        ...

    def is_isomorphic_vf2(self, graph: Graph) -> bool:
        """Compare isomorphism by VF2 algorithm,
        one of the high performance isomorphic algorithms.
        """
        ...

    def is_isomorphic_degree_code(self, graph: Graph) -> bool:
        """Compare isomorphism by degree code algorithm.

        + <https://doi.org/10.1115/1.2919236>
        """
        ...

    def duplicate(self, vertices: Iterable[int], times: int) -> Graph:
        """Make graph duplicate by specific `vertices`. Return a new graph."""
        ...

    def copy(self) -> Graph:
        """The copy method of the Graph object."""
        ...
