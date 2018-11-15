# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class.

The algorithms references:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
import sys
from typing import (
    Sequence,
    Tuple,
    Dict,
    Iterator,
)
from itertools import combinations
from libcpp.map cimport map


@cython.final
cdef class Graph:

    """NetworkX-like graph class."""

    def __cinit__(self, edges: Sequence[Tuple[int, int]]):
        """edges: [(l1, l2), ...]"""
        self.edges = tuple(edges)

        # nodes
        cdef int p1, p2
        cdef list nodes = []
        for p1, p2 in self.edges:
            if p1 not in nodes:
                nodes.append(p1)
            if p2 not in nodes:
                nodes.append(p2)
        self.nodes = tuple(nodes)

        # adj
        cdef int n
        self.adj = {n: self.neighbors(n) for n in self.nodes}

    cdef inline tuple neighbors(self, int n):
        """Neighbors except the node."""
        cdef list neighbors = []
        cdef int l1, l2
        for l1, l2 in self.edges:
            if n == l1:
                neighbors.append(l2)
            if n == l2:
                neighbors.append(l1)
        return tuple(neighbors)

    cpdef list degrees(self):
        """Return number of neighbors par node."""
        cdef int n
        return [(n, len(neighbors)) for n, neighbors in self.adj.items()]

    cpdef int dof(self):
        """Return degrees of freedom."""
        return 3 * (len(self.nodes) - 1) - (2 * len(self.edges))

    cpdef bool is_connected(self, int with_out = -1):
        """Return True if the graph is not isolated."""
        cdef int neighbors
        cdef int index = 0
        cdef list nodes = []
        # Change start node if index zero has been excluded.
        if with_out == self.nodes[0]:
            nodes.append(self.nodes[1])
        else:
            nodes.append(self.nodes[0])

        # Search node by node.
        while index < len(nodes):
            for neighbors in self.adj[nodes[index]]:
                if (neighbors not in nodes) and (neighbors != with_out):
                    nodes.append(neighbors)
            index += 1

        if with_out != -1:
            nodes.append(with_out)
        return len(nodes) == len(self.nodes)

    cpdef bool has_cut_link(self):
        """Return True if the graph has any cut links."""
        cdef int node, degree
        for node, degree in self.degrees():
            # Only for multiple links.
            if degree > 2:
                # Remove a multiple link should be keep connected.
                if not self.is_connected(node):
                    return True
        return False

    cpdef bool is_degenerate(self):
        """Return True if this kinematic chain is degenerate.
        
        + Prue all multiple contracted links recursively.
        + Check the DOF of sub-graph if it is lower then zero.
        """
        if self.has_triangles():
            return True
        cdef int n1, n2
        cdef tuple neighbors
        cdef list c_l
        cdef Graph g = self.copy()
        cdef set mcl = set()
        while True:
            mcl.clear()
            mcl.update(g.multi_contracted_links())
            if not mcl:
                break
            c_l = []
            for n1, n2 in g.edges:
                if not {n1, n2} & mcl:
                    c_l.append((n1, n2))
            g = Graph(c_l)
            if {len(neighbors) for neighbors in g.adj.values()} == {2}:
                # The graph is a basic loop.
                break

        return g.dof() < 1

    cpdef bool is_isomorphic(self, Graph graph):
        """Return True if two graphs is isomorphic."""
        cdef GraphMatcher gm_gh = GraphMatcher(self, graph)
        return gm_gh.is_isomorphic()

    cpdef bool is_planar(self):
        """Return True if the graph is planar."""
        # TODO: check_planarity.
        return True

    cdef list link_types(self):
        """Return types of each link."""
        cdef int n, d
        return sorted([d for n, d in self.degrees()])

    cdef int node_distance(self, int u, int v):
        """Find the distance between u and v."""
        if v in self.adj[u]:
            return 1
        return 0

    cdef list multi_contracted_links(self):
        """Return a list of multiple contracted links."""
        cdef int n1, n2, index, neighbor
        cdef list m_links
        for n1 in self.nodes:
            # Only for binary link.
            if len(self.adj[n1]) != 2:
                continue
            # n1 is not collected yet.
            index = 0
            m_links = [n1]
            while index < len(m_links):
                for neighbor in self.adj[m_links[index]]:
                    if len(self.adj[neighbor]) == 2:
                        if neighbor not in m_links:
                            m_links.append(neighbor)
                index += 1
            if len(m_links) > 1:
                return m_links
        return []

    cdef bool has_triangles(self):
        """Return True if the graph has triangles."""
        cdef int n1, n2
        cdef tuple neighbors
        for neighbors in self.adj.values():
            for n1 in neighbors:
                for n2 in neighbors:
                    if n1 == n2:
                        continue
                    if n1 in self.adj[n2]:
                        return True
        return False

    cpdef Graph copy(self):
        """Copy self."""
        return Graph(self.edges)

    cpdef Graph subgraph(self, tuple nodes):
        """Return a sub-graph from self."""
        cdef set nodes_set = set(nodes)
        cdef tuple edge
        return Graph([edge for edge in self.edges if set(edge) <= nodes_set])


@cython.final
cdef class GraphMatcher:

    """GraphMatcher and GMState class from NetworkX.

    Copyright (C) 2007-2009 by the NetworkX maintainers
    All rights reserved.
    BSD license.

    This work was originally coded by Christopher Ellison
    as part of the Computational Mechanics Python (CMPy) project.
    James P. Crutchfield, principal investigator.
    Complexity Sciences Center and Physics Department, UC Davis.
    """

    cdef Graph g1, g2
    cdef set g1_nodes, g2_nodes
    cdef dict core_1, core_2, inout_1, inout_2, mapping
    cdef GMState state

    def __cinit__(self, g1: Graph, g2: Graph):
        self.g1 = g1
        self.g2 = g2
        self.g1_nodes = set(g1.nodes)
        self.g2_nodes = set(g2.nodes)

        # Set recursion limit.
        cdef int old_recursion_limit = sys.getrecursionlimit()
        cdef int expected_max_recursion_level = len(self.g2.nodes)
        if old_recursion_limit < 1.5 * expected_max_recursion_level:
            # Give some breathing room.
            sys.setrecursionlimit(int(1.5 * expected_max_recursion_level))

        # Initialize state
        self.initialize()

    cdef inline void initialize(self):
        """Re-initializes the state of the algorithm."""
        # core_1[n] contains the index of the node paired with n, which is m,
        #            provided n is in the mapping.
        # core_2[m] contains the index of the node paired with m, which is n,
        #            provided m is in the mapping.
        self.core_1 = {}
        self.core_2 = {}

        # See the paper for definitions of M_x and T_x^{y}
        # inout_1[n]  is non-zero if n is in M_1 or in T_1^{inout}
        # inout_2[m]  is non-zero if m is in M_2 or in T_2^{inout}
        #
        # The value stored is the depth of the SSR tree when the node became
        # part of the corresponding set.
        self.inout_1 = {}
        self.inout_2 = {}
        # Practically, these sets simply store the nodes in the subgraph.

        self.state = GMState(self)

        # Provide a convenient way to access the isomorphism mapping.
        self.mapping = self.core_1.copy()

    def candidate_pairs_iter(self) -> Iterator[Tuple[int, int]]:
        """Iterator over candidate pairs of nodes in g1 and g2."""
        cdef int node
        # First we compute the inout-terminal sets.
        cdef set s1 = set(self.inout_1) - set(self.core_1)
        cdef set s2 = set(self.inout_2) - set(self.core_2)
        cdef list t1_inout = [node for node in self.g1_nodes if (node in s1)]
        cdef list t2_inout = [node for node in self.g2_nodes if (node in s2)]
        # If t1_inout and t2_inout are both nonempty.
        # P(s) = t1_inout x {min t2_inout}
        if t1_inout and t2_inout:
            for node in t1_inout:
                yield node, min(t2_inout)
        else:
            # If t1_inout and t2_inout were both empty....
            # P(s) = (N_1 - M_1) x {min (N_2 - M_2)}
            # # if not (t1_inout or t2_inout):
            # as suggested by  [2], incorrect
            # as inferred from [1], correct
            # First we determine the candidate node for g2
            for node in self.g1.nodes:
                if node not in self.core_1:
                    yield node, min(self.g2_nodes - set(self.core_2))
        # For all other cases, we don't have any candidate pairs.

    cdef bool is_isomorphic(self):
        """Returns True if g1 and g2 are isomorphic graphs."""
        # Let's do two very quick checks!
        # QUESTION: Should we call faster_graph_could_be_isomorphic(g1,g2)?
        # For now, graph3 just copy the code.

        # Check global properties
        if len(self.g1.nodes) != len(self.g2.nodes):
            return False

        # Check local properties
        if self.g1.link_types() != self.g2.link_types():
            return False
        try:
            next(self.isomorphisms_iter())
            return True
        except StopIteration:
            return False

    def isomorphisms_iter(self) -> Iterator[Dict[int, Tuple[int]]]:
        """Generator over isomorphisms between g1 and g2.

        Declare that we are looking for a graph-graph isomorphism.
        """
        self.initialize()
        yield from self.match()

    def match(self) -> Iterator[Dict[int, Tuple[int]]]:
        """Extends the isomorphism mapping."""
        cdef int g1_node, g2_node
        cdef GMState new_state
        if len(self.core_1) == len(self.g2.nodes):
            # Save the final mapping, otherwise garbage collection deletes it.
            self.mapping = self.core_1.copy()
            # The mapping is complete.
            yield self.mapping
        else:
            for g1_node, g2_node in self.candidate_pairs_iter():
                if self.syntactic_feasibility(g1_node, g2_node):
                    # Recursive call, adding the feasible state.
                    new_state = self.state.__class__(self, g1_node, g2_node)

                    yield from self.match()

                    # restore data structures
                    new_state.restore()

    cdef inline bool syntactic_feasibility(self, int g1_node, int g2_node):
        """Returns True if adding (g1_node, g2_node) is syntactically feasible.
        
        The VF2 algorithm was designed to work with graphs having, at most,
        one edge connecting any two nodes.  This is not the case when
        dealing with an MultiGraphs.
        
        Basically, when we test the look-ahead rules R_neighbor, we will
        make sure that the number of edges are checked. We also add
        a R_self check to verify that the number of self loops is acceptable.
        
        Users might be comparing Graph instances with MultiGraph instances.
        So the generic GraphMatcher class must work with MultiGraphs.
        Care must be taken since the value in the innermost dictionary is a
        singlet for Graph instances.  For MultiGraphs, the value in the
        innermost dictionary is a list.
        """

        # ## Test at each step to get a return value as soon as possible.

        # ## Look ahead 0
        # R_self
        # The number of self loops for g1_node must equal the number of
        # self-loops for g2_node. Without this check, we would fail on
        # R_neighbor at the next recursion level. But it is good to prune the
        # search tree now.
        if self.g1.node_distance(g1_node, g1_node) != self.g2.node_distance(g2_node, g2_node):
            return False
        # R_neighbor
        # For each neighbor n' of n in the partial mapping, the corresponding
        # node m' is a neighbor of m, and vice versa. Also, the number of
        # edges must be equal.
        cdef int neighbor
        for neighbor in self.g1.adj[g1_node]:
            if neighbor in self.core_1:
                if self.core_1[neighbor] not in self.g2.adj[g2_node]:
                    return False
                elif (
                        self.g1.node_distance(neighbor, g1_node) !=
                        self.g2.node_distance(self.core_1[neighbor], g2_node)
                ):
                    return False
        for neighbor in self.g2.adj[g2_node]:
            if neighbor in self.core_2:
                if self.core_2[neighbor] not in self.g1.adj[g1_node]:
                    return False
                elif (
                        self.g1.node_distance(self.core_2[neighbor], g1_node) !=
                        self.g2.node_distance(neighbor, g2_node)
                ):
                    return False

        # ## Look ahead 1
        # R_term_inout
        # The number of neighbors of n that are in T_1^{inout} is equal to the
        # number of neighbors of m that are in T_2^{inout}, and vice versa.
        cdef int num1 = 0
        for neighbor in self.g1.adj[g1_node]:
            if (neighbor in self.inout_1) and (neighbor not in self.core_1):
                num1 += 1
        cdef int num2 = 0
        for neighbor in self.g2.adj[g2_node]:
            if (neighbor in self.inout_2) and (neighbor not in self.core_2):
                num2 += 1
        if num1 != num2:
            return False

        # ## Look ahead 2
        # R_new
        # The number of neighbors of n that are neither in the core_1 nor
        # T_1^{inout} is equal to the number of neighbors of m
        # that are neither in core_2 nor T_2^{inout}.
        num1 = 0
        for neighbor in self.g1.adj[g1_node]:
            if neighbor not in self.inout_1:
                num1 += 1
        num2 = 0
        for neighbor in self.g2.adj[g2_node]:
            if neighbor not in self.inout_2:
                num2 += 1
        return num1 == num2


@cython.final
cdef class GMState:

    cdef GraphMatcher gm
    cdef int g1_node, g2_node, depth

    def __cinit__(
        self,
        gm: GraphMatcher,
        g1_node: int = -1,
        g2_node: int = -1
    ):
        """Initializes GMState object.

        Pass in the GraphMatcher to which this GMState belongs and the
        new node pair that will be added to the GraphMatcher's current
        isomorphism mapping.
        """
        self.gm = gm

        # Initialize the last stored node pair.
        self.g1_node = -1
        self.g2_node = -1
        self.depth = len(gm.core_1)

        if g1_node == -1 or g2_node == -1:
            # Then we reset the class variables
            gm.core_1 = {}
            gm.core_2 = {}
            gm.inout_1 = {}
            gm.inout_2 = {}
        cdef set new_nodes
        cdef int node, neighbor
        # Watch out! g1_node == 0 should evaluate to True.
        if g1_node != -1 and g2_node != -1:
            # Add the node pair to the isomorphism mapping.
            gm.core_1[g1_node] = g2_node
            gm.core_2[g2_node] = g1_node

            # Store the node that was added last.
            self.g1_node = g1_node
            self.g2_node = g2_node

            # Now we must update the other two vectors.
            # We will add only if it is not in there already!
            self.depth = len(gm.core_1)

            # First we add the new nodes...
            if g1_node not in gm.inout_1:
                gm.inout_1[g1_node] = self.depth
            if g2_node not in gm.inout_2:
                    gm.inout_2[g2_node] = self.depth

            # Now we add every other node...

            # Updates for T_1^{inout}
            new_nodes = set()
            for node in gm.core_1:
                new_nodes.update([neighbor for neighbor in gm.g1.adj[node] if neighbor not in gm.core_1])
            for node in new_nodes:
                if node not in gm.inout_1:
                    gm.inout_1[node] = self.depth

            # Updates for T_2^{inout}
            new_nodes = set()
            for node in gm.core_2:
                new_nodes.update([neighbor for neighbor in gm.g2.adj[node] if neighbor not in gm.core_2])
            for node in new_nodes:
                if node not in gm.inout_2:
                    gm.inout_2[node] = self.depth

    cdef void restore(self):
        """Deletes the GMState object and restores the class variables."""
        # First we remove the node that was added from the core vectors.
        # Watch out! g1_node == 0 should evaluate to True.
        if self.g1_node != -1 and self.g2_node != -1:
            del self.gm.core_1[self.g1_node]
            del self.gm.core_2[self.g2_node]

        # Now we revert the other two vectors.
        # Thus, we delete all entries which have this depth level.
        cdef dict vector
        cdef int node
        for vector in (self.gm.inout_1, self.gm.inout_2):
            for node in tuple(vector):
                if vector[node] == self.depth:
                    del vector[node]
