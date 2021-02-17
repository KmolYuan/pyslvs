# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Graph class.

The algorithm reference:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from libcpp.pair cimport pair
import sys
from typing import Tuple, Dict, Iterator
from itertools import permutations, groupby
from numpy import zeros, uint8 as np_uint

ctypedef pair[int, int] ipair


cpdef list link_assortment(Graph g):
    """Return link assortment of the graph."""
    assortment = [0]
    if not g.edges:
        return assortment
    cdef int d, n
    for n in g.degrees().values():
        if n < 2:
            continue
        d = n - 2
        while d >= len(assortment):
            assortment.append(0)
        assortment[d] += 1
    return assortment


cpdef list contracted_link_assortment(Graph g):
    """Return contracted link assortment of the graph."""
    if not g.edges:
        return [0]
    assortment = [0] * link_assortment(g)[0]
    counted = set()
    for mcl in _multi_contracted_links(g, False):
        counted.update(mcl)
        assortment[len(mcl) - 1] += 1
    # For single contracted links
    cdef int n, d
    for n, d in g.degrees().items():
        if d != 2:
            continue
        if n not in counted:
            assortment[0] += 1
    return assortment


cdef list _multi_contracted_links(Graph g, bint only_one):
    """Return a list of multiple contracted links."""
    cdef int n1, n2, index, neighbor
    contracted_links = []
    counted = set()
    for n1 in g.vertices:
        # Only for binary link.
        if len(g.adj[n1]) != 2:
            continue
        # n1 is not collected yet.
        if n1 in counted:
            continue
        index = 0
        c_links = [n1]
        while index < len(c_links):
            for neighbor in g.adj[c_links[index]]:
                if len(g.adj[neighbor]) == 2 and neighbor not in c_links:
                    c_links.append(neighbor)
            index += 1
        if len(c_links) > 1:
            if only_one:
                return c_links
            counted.update(c_links)
            contracted_links.append(tuple(c_links))
    if only_one:
        return []
    else:
        return contracted_links


cpdef list labeled_enumerate(Graph g):
    """Enumerate each node with labeled except isomorphism."""
    result = []
    cdef int n1, n2
    cdef Graph g1, g2
    for n1 in g.vertices:
        g1 = Graph.__new__(Graph, [e for e in g.edges if n1 not in e])
        for n2, g2 in result:
            if g1.is_isomorphic(g2):
                break
        else:
            result.append((n1, g1))
    return result


@cython.final
cdef class Graph:
    """The undirected graph class, support multigraph."""

    def __cinit__(self, object edges):
        self.edges = tuple(edges)
        # vertices
        cdef int p1, p2
        vertices = []
        for p1, p2 in self.edges:
            if p1 not in vertices:
                vertices.append(p1)
            if p2 not in vertices:
                vertices.append(p2)
        self.vertices = tuple(vertices)
        # adj
        self.adj = {n: self.neighbors(n) for n in self.vertices}

    cpdef void add_vertices(self, object vertices):
        """Add vertices from iterable object `vertices`."""
        self.vertices = tuple(set(self.vertices) | set(vertices))

    cpdef void add_edge(self, int n1, int n2):
        """Add edge `n1` to `n2`."""
        self.edges += ((n1, n2),)
        self.vertices = tuple(set(self.vertices) | {n1, n2})
        cdef int n
        for n in (n1, n2):
            self.adj[n] = self.neighbors(n)

    cpdef void add_path(self, object new_nodes):
        edges = list(self.edges)
        vertices = set(self.vertices)
        cdef int n1 = -1
        cdef int n2
        for n2 in new_nodes:
            if n1 != -1:
                edges.append((n1, n2))
                vertices.update((n1, n2))
            n1 = n2
        self.edges = tuple(edges)
        self.vertices = tuple(vertices)
        for n1 in vertices:
            self.adj[n1] = self.neighbors(n1)

    cpdef void remove_edge(self, int n1, int n2):
        """Remove edge(s) {n1, n2} once if exist, otherwise do nothing."""
        vertices = set()
        edges = []
        cdef bint once = False
        for edge in self.edges:
            if {n1, n2} != set(edge) or once:
                vertices.update(edge)
                edges.append(edge)
            else:
                once = True
        self.edges = tuple(edges)
        self.vertices = tuple(vertices)
        cdef int n
        for n in self.vertices:
            edge = self.neighbors(n)
            if edge:
                self.adj[n] = edge
            else:
                self.adj.pop(n, None)

    cpdef int dof(self):
        """Return DOF of the graph.

        !!! note
            DOF is the Degree of Freedoms to a mechanism.

            In the [Graph] objects, all vertices will assumed as revolute
            joints (1 DOF).

            $$
            F = 3(N_L - 1) - 2N_J
            $$
        """
        return 3 * (len(self.vertices) - 1) - (2 * len(self.edges))

    cpdef inline tuple neighbors(self, int n):
        """Return the neighbors of the vertex `n`."""
        neighbors = []
        cdef int l1, l2
        for l1, l2 in self.edges:
            if n == l1:
                neighbors.append(l2)
            if n == l2:
                neighbors.append(l1)
        return tuple(neighbors)

    cpdef dict degrees(self):
        """Return the degrees of each vertex."""
        return {n: len(neighbors) for n, neighbors in self.adj.items()}

    cpdef ullong degree_code(self):
        """Generate a degree code.

        With a sorted vertices mapping by the degrees of each vertex,
        regenerate a new adjacency matrix.
        A binary code can be found by concatenating the upper right elements.
        The degree code is the maximum value of the permutation.
        """
        if len(self.vertices) < 2:
            return 0
        # Create a new mapping
        degrees = self.degrees()
        # Permute the group to find the max degree code
        # "per1" is the mapping of the graph
        # "prefix" is the last choose of the candidate (can be multiple)
        per1 = []
        prefix = []
        cdef int i, n1, n2
        cdef ullong code, sub_code
        for _, g in groupby(
            sorted(degrees, key=degrees.__getitem__, reverse=True),
            key=degrees.__getitem__
        ):
            code = 0
            order = None
            g = tuple(g)
            for pre in iter(prefix if len(prefix) > 1 else [()]):
                for per2 in permutations(g):
                    if len(per2) == 1:
                        order = pre + per2
                        prefix = [order]
                        break
                    sub_code = 0
                    # Calculate sub code
                    per3 = tuple(per1) + pre + per2
                    for i, n1 in enumerate(per3):
                        for n2 in per3[i + 1:]:
                            sub_code <<= 1
                            sub_code += n2 in self.adj[n1]
                    # Compare sub code
                    if sub_code > code or order is None:
                        code = sub_code
                        order = pre + per2
                        prefix = [order]
                    elif sub_code == code:
                        # Will affect the subsequent arrangement
                        prefix.append(pre + per2)
            if len(prefix) == 1:
                per1.extend(order)
        # Check the last one prefix candidate
        if not set(order) <= set(per1):
            per1.extend(order)
        # Calculate the degree code
        code = 0
        for i, n1 in enumerate(per1):
            for n2 in per1[i + 1:]:
                code <<= 1
                code += n2 in self.adj[n1]
        return code

    cpdef double[:, :] adjacency_matrix(self):
        """Generate a adjacency matrix.

        Assume the matrix $A[i, j] = A[j, i]$.
        Where $A[i, j] = 1$ if edge `(i, j)` exist.
        """
        cdef int n = len(self.vertices)
        cdef double[:, :] am = zeros((n, n), dtype=np_uint)
        cdef int n1, n2
        for n1, n2 in self.edges:
            am[n1, n2] += 1
            am[n2, n1] += 1
        return am

    cpdef bint is_connected(self, int without = -1):
        """Return `True` if the graph is connected.
        Set the argument `without` to ignore one vertex.
        """
        if not self.vertices:
            return True
        cdef int neighbors
        cdef int index = 0
        vertices = []
        # Change start node if index zero has been excluded.
        if without == self.vertices[0]:
            vertices.append(self.vertices[1])
        else:
            vertices.append(self.vertices[0])
        # Search node by node.
        while index < len(vertices):
            for neighbors in self.adj[vertices[index]]:
                if (neighbors not in vertices) and neighbors != without:
                    vertices.append(neighbors)
            index += 1
        if without != -1:
            vertices.append(without)
        return len(vertices) == len(self.vertices)

    cpdef bint has_cut_link(self):
        """Return true if the graph has any cut links."""
        cdef int n, d
        for n, d in self.degrees().items():
            # Only for multiple links.
            if d > 2:
                # Remove a multiple link should be keep connected.
                if not self.is_connected(n):
                    return True
        return False

    cpdef bint is_degenerate(self):
        """Return true if this kinematic chain is degenerate.

        + Prue all multiple contracted links recursively.
        + Check the DOF of sub-graph if it is lower then zero.
        """
        if self.has_triangle():
            return True
        cdef int n1, n2
        cdef Graph g = self.copy()
        mcl = set()
        while True:
            mcl.update(_multi_contracted_links(g, True))
            if not mcl:
                break
            c_l = []
            for n1, n2 in g.edges:
                if not {n1, n2} & mcl:
                    c_l.append((n1, n2))
            mcl.clear()
            # Pruned graph
            g = Graph.__new__(Graph, c_l)
            if {n2 for n2 in g.degrees().values()} == {2}:
                # The graph is a basic loop.
                break
        # Check the DOF
        return g.dof() < 1

    cpdef bint has_triangle(self):
        """Return true if the graph has triangle."""
        cdef int n1, n2
        for neighbors in self.adj.values():
            for n1 in neighbors:
                for n2 in neighbors:
                    if n1 == n2:
                        continue
                    if n1 in self.adj[n2]:
                        return True
        return False

    cpdef bint is_isomorphic(self, Graph g):
        """Return true if two graphs is isomorphic.

        Default is using VF2 algorithm.
        """
        return self.is_isomorphic_vf2(g)

    cpdef bint is_isomorphic_vf2(self, Graph g):
        """Compare isomorphism by VF2 algorithm,
        one of the high performance isomorphic algorithms.
        """
        cdef GraphMatcher gm_gh = GraphMatcher.__new__(GraphMatcher, self, g)
        return gm_gh.is_isomorphic()

    cpdef bint is_isomorphic_degree_code(self, Graph g):
        """Compare isomorphism by degree code algorithm.

        + <https://doi.org/10.1115/1.2919236>
        """
        return self.degree_code() == g.degree_code()

    cpdef Graph duplicate(self, object vertices, int times):
        """Make graph duplicate by specific `vertices`. Return a new graph."""
        if times < 1:
            raise ValueError("please input a number larger than 1.")
        cdef int max_num = max(self.vertices) + 1
        mapping = {}
        cdef int i, n1, n2
        edges = set(self.edges)
        for i in range(times):
            for n1 in sorted(set(vertices)):
                mapping[n1] = max_num
                max_num += 1
            for n1, n2 in self.edges:
                if n1 in mapping:
                    n1 = mapping[n1]
                if n2 in mapping:
                    n2 = mapping[n2]
                if n1 > n2:
                    n1, n2 = n2, n1
                edges.add((n1, n2))
        return Graph.__new__(Graph, edges)

    cpdef Graph copy(self):
        """The copy method of the Graph object."""
        return Graph.__new__(Graph, self.edges)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({list(self.edges)})"


cdef bint _is_adjacent(Graph g, int u, int v):
    """Find the distance between u and v."""
    if v in g.adj[u]:
        return True
    return False


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

    def __cinit__(self, Graph g1, Graph g2):
        self.g1 = g1
        self.g2 = g2
        self.g1_nodes = set(g1.vertices)
        self.g2_nodes = set(g2.vertices)
        # Set recursion limit
        cdef int old_recursion_limit = sys.getrecursionlimit()
        cdef int expected_max_recursion_level = len(self.g2.vertices)
        if old_recursion_limit < 1.5 * expected_max_recursion_level:
            # Give some breathing room
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
        self.state = GMState.__new__(GMState, self)
        # Provide a convenient way to access the isomorphism mapping.
        self.mapping = self.core_1.copy()

    def candidate_pairs_iter(self) -> Iterator[Tuple[int, int]]:
        """Iterator over candidate pairs of nodes in g1 and g2."""
        # First we compute the inout-terminal sets.
        s1 = set(self.inout_1) - set(self.core_1)
        s2 = set(self.inout_2) - set(self.core_2)
        t1_inout = [node for node in self.g1_nodes if node in s1]
        t2_inout = [node for node in self.g2_nodes if node in s2]
        # If t1_inout and t2_inout are both nonempty.
        # P(s) = t1_inout x {min t2_inout}
        cdef int node
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
            for node in self.g1.vertices:
                if node not in self.core_1:
                    yield node, min(self.g2_nodes - set(self.core_2))
        # For all other cases, we don't have any candidate pairs.

    cdef bint is_isomorphic(self):
        """Returns True if g1 and g2 are isomorphic graphs."""
        # Let's do two very quick checks!
        # QUESTION: Should we call faster_graph_could_be_isomorphic(g1,g2)?
        # For now, graph3 just copy the code.
        # Check global properties
        if len(self.g1.vertices) != len(self.g2.vertices):
            return False
        # Check local properties
        if (
            sorted([len(neighbor) for neighbor in self.g1.adj.values()])
            != sorted([len(neighbor) for neighbor in self.g2.adj.values()])
        ):
            return False
        self.initialize()
        try:
            next(self.match())
            return True
        except StopIteration:
            return False

    def match(self) -> Iterator[Dict[int, Tuple[int]]]:
        """Extends the isomorphism mapping."""
        cdef int g1_node, g2_node
        cdef GMState new_state
        if len(self.core_1) == len(self.g2.vertices):
            # Save the final mapping, otherwise garbage collection deletes it.
            self.mapping = self.core_1.copy()
            # The mapping is complete.
            yield self.mapping
        else:
            for g1_node, g2_node in self.candidate_pairs_iter():
                if self.syntactic_feasibility(g1_node, g2_node):
                    # Recursive call, adding the feasible state.
                    new_state = type(self.state)(self, g1_node, g2_node)
                    yield from self.match()
                    # restore data structures
                    new_state.restore()

    cdef inline bint syntactic_feasibility(self, int g1_node, int g2_node):
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
        if (
            _is_adjacent(self.g1, g1_node, g1_node)
            != _is_adjacent(self.g2, g2_node, g2_node)
        ):
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
                        _is_adjacent(self.g1, neighbor, g1_node)
                        != _is_adjacent(self.g2, self.core_1[neighbor], g2_node)
                ):
                    return False
        for neighbor in self.g2.adj[g2_node]:
            if neighbor in self.core_2:
                if self.core_2[neighbor] not in self.g1.adj[g1_node]:
                    return False
                elif (
                        _is_adjacent(self.g1, self.core_2[neighbor], g1_node)
                        != _is_adjacent(self.g2, neighbor, g2_node)
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
    """Graph matcher state."""
    cdef GraphMatcher gm
    cdef int g1_node, g2_node, depth

    def __cinit__(
        self,
        GraphMatcher gm,
        int g1_node=-1,
        int g2_node=-1
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

        cdef int node
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
                new_nodes.update([
                    neighbor for neighbor in gm.g1.adj[node]
                    if neighbor not in gm.core_1
                ])
            for node in new_nodes:
                if node not in gm.inout_1:
                    gm.inout_1[node] = self.depth

            # Updates for T_2^{inout}
            new_nodes = set()
            for node in gm.core_2:
                new_nodes.update([
                    neighbor for neighbor in gm.g2.adj[node]
                    if neighbor not in gm.core_2
                ])
            for node in new_nodes:
                if node not in gm.inout_2:
                    gm.inout_2[node] = self.depth

    cdef void restore(self):
        """Deletes the GMState object and restores the class variables."""
        # First we remove the node that was added from the core vectors.
        # Watch out! g1_node == 0 should evaluate to True.
        if self.g1_node != -1 and self.g2_node != -1:
            self.gm.core_1.pop(self.g1_node)
            self.gm.core_2.pop(self.g2_node)

        # Now we revert the other two vectors.
        # Thus, we delete all entries which have this depth level.
        cdef int node
        for vector in (self.gm.inout_1, self.gm.inout_2):
            for node in tuple(vector):
                if vector[node] == self.depth:
                    vector.pop(node)
