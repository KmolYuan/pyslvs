# -*- coding: utf-8 -*-
# cython: language_level=3

"""Type synthesis."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

cimport cython
from typing import (
    Sequence,
    Tuple,
    Dict,
    Iterator,
)
from itertools import combinations, product
import sys
from time import time
from cpython cimport bool
from numpy cimport ndarray
from numpy import (
    zeros as np_zeros,
    int32 as np_int32,
)


@cython.final
cdef class Graph:
    
    """NetworkX-like graph class."""
    
    cdef public tuple edges
    cdef tuple nodes
    cdef dict adj
    
    def __cinit__(self, edges: Sequence[Tuple[int, int]]):
        # edges
        """edges: ((l1, l2), ...)"""
        self.edges = tuple(edges)
        # nodes
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
    
    cpdef Graph compose(self, Graph graph):
        return Graph(set(self.edges) | set(graph.edges))
    
    cpdef bool out_of_limit(self, ndarray limit):
        cdef int n
        for n in self.adj:
            if len(self.adj[n]) > limit[n]:
                return True
        return False
    
    cpdef bool has_triangles(self):
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
    
    cpdef bool is_connected(self):
        cdef int index = 0
        cdef list nodes = [self.nodes[index]]
        while index < len(nodes):
            for neighbor in self.adj[nodes[index]]:
                if neighbor not in nodes:
                    nodes.append(neighbor)
            index += 1
        return len(nodes) == len(self.nodes)
    
    cpdef bool is_isomorphic(self, Graph graph):
        cdef GraphMatcher gm_gh = GraphMatcher(self, graph)
        return gm_gh.is_isomorphic()
    
    cpdef list links(self):
        return sorted([len(neighbors) for neighbors in self.adj.values()])
    
    cpdef int number_of_edges(self, int u, int v):
        if v in self.adj[u]:
            return 1
        return 0


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
    
    # Re-initializes the state of the algorithm.
    cdef inline void initialize(self):
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
    
    # Generator candidate_pairs_iter()
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
    
    # Returns True if g1 and g2 are isomorphic graphs.
    cpdef bool is_isomorphic(self):
        # Let's do two very quick checks!
        # QUESTION: Should we call faster_graph_could_be_isomorphic(g1,g2)?
        # For now, graph3 just copy the code.
        
        # Check global properties
        if len(self.g1.nodes) != len(self.g2.nodes):
            return False
        
        # Check local properties
        if self.g1.links() != self.g2.links():
            return False
        try:
            next(self.isomorphisms_iter())
            return True
        except StopIteration:
            return False
    
    # Generator isomorphisms_iter()
    # Generator over isomorphisms between g1 and g2.
    def isomorphisms_iter(self) -> Iterator[Dict[int, Tuple[int]]]:
        # Declare that we are looking for a graph-graph isomorphism.
        self.initialize()
        cdef dict mapping
        for mapping in self.match():
            yield mapping
    
    # Generator match()
    # Extends the isomorphism mapping.
    def match(self) -> Iterator[Dict[int, Tuple[int]]]:
        cdef int g1_node, g2_node
        cdef GMState new_state
        cdef dict mapping
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
                    for mapping in self.match():
                        yield mapping
                    # restore data structures
                    new_state.restore()
    
    # Returns True if adding (g1_node, g2_node) is syntactically feasible.
    cdef inline bool syntactic_feasibility(self, int g1_node, int g2_node):
        # The VF2 algorithm was designed to work with graphs having, at most,
        # one edge connecting any two nodes.  This is not the case when
        # dealing with an MultiGraphs.
        # 
        # Basically, when we test the look-ahead rules R_neighbor, we will
        # make sure that the number of edges are checked. We also add
        # a R_self check to verify that the number of self loops is acceptable.
        # 
        # Users might be comparing Graph instances with MultiGraph instances.
        # So the generic GraphMatcher class must work with MultiGraphs.
        # Care must be taken since the value in the innermost dictionary is a
        # singlet for Graph instances.  For MultiGraphs, the value in the
        # innermost dictionary is a list.
        
        # ## Test at each step to get a return value as soon as possible.
        
        # ## Look ahead 0
        # R_self
        # The number of self loops for g1_node must equal the number of
        # self-loops for g2_node. Without this check, we would fail on
        # R_neighbor at the next recursion level. But it is good to prune the
        # search tree now.
        if self.g1.number_of_edges(g1_node, g1_node) != self.g2.number_of_edges(g2_node, g2_node):
            return False
        # R_neighbor
        # For each neighbor n' of n in the partial mapping, the corresponding
        # node m' is a neighbor of m, and vice versa. Also, the number of
        # edges must be equal.
        cdef int neighbor
        for neighbor in self.g1.adj[g1_node]:
            if neighbor in self.core_1:
                if not (self.core_1[neighbor] in self.g2.adj[g2_node]):
                    return False
                elif self.g1.number_of_edges(neighbor, g1_node) != self.g2.number_of_edges(self.core_1[neighbor], g2_node):
                    return False
        for neighbor in self.g2.adj[g2_node]:
            if neighbor in self.core_2:
                if not (self.core_2[neighbor] in self.g1.adj[g1_node]):
                    return False
                elif self.g1.number_of_edges(self.core_2[neighbor], g1_node) != self.g2.number_of_edges(neighbor, g2_node):
                    return False
        
        # ## Look ahead 1
        # R_terminout
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
    
    cpdef void restore(self):
        """Deletes the GMState object and restores the class variables."""
        # First we remove the node that was added from the core vectors.
        # Watch out! g1_node == 0 should evaluate to True.
        if self.g1_node!=-1 and self.g2_node!=-1:
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


cdef inline bool is_isomorphic(Graph graph1, list answer):
    """Return True if the graph is isomorphic with list."""
    cdef Graph graph2
    for i, graph2 in enumerate(answer):
        if graph1.is_isomorphic(graph2):
            return True
    return False


cdef inline list connection_get(int i, tuple connection):
    cdef tuple c
    return [c for c in connection if (i in c)]


cpdef tuple topo(
    object link_num,
    bool degenerate = True,
    object job_func = None,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    link_num = [L2, L3, L4, ...]
    links = [[number_code]: joint_number, ...]
    """
    # Initial time.
    cdef double t0 = time()
    
    # Number of all joints.
    cdef int joint_count = sum(link_num)
    
    # Number of joins in each links.
    cdef ndarray links = np_zeros((joint_count,), dtype=np_int32)
    cdef int i, j, t, name
    for i in range(joint_count):
        name = i
        for j, t in enumerate(link_num):
            if i < t:
                links[name] = j + 2
                break
            i -= t
        else:
            raise RuntimeError("AAA")
    
    # connections = [(1, 2), (1, 3), ..., (2, 3), (2, 4), ...]
    cdef tuple connections = tuple(combinations(range(joint_count), 2))
    
    # Result list.
    cdef list edges_combinations = []
    cdef list matched = []
    
    # Find ALL results.
    cdef int link, count, progress_value
    cdef list match
    cdef tuple m
    cdef Graph graph1, graph2, graph3
    for link, count in enumerate(links):
        # Other of joints that the link connect with.
        match = [Graph(m) for m in combinations(
            connection_get(link, connections),
            count
        )]
        
        if not edges_combinations:
            edges_combinations = match
            continue
        
        if job_func:
            progress_value = len(edges_combinations) * len(match)
            job_func(
                f"Match link # {link} / {len(links) - 1}\n"
                f"Possibilities: {progress_value}",
                progress_value
            )
        
        matched.clear()
        for graph1, graph2 in product(edges_combinations, match):
            if stop_func and stop_func():
                break
            graph3 = graph1.compose(graph2)
            # Out of limit.
            if graph3.out_of_limit(links):
                continue
            # Has triangles.
            if degenerate and graph3.has_triangles():
                continue
            # TODO: Check isomorphic here?
            #if is_isomorphic(graph3, matched):
            #    continue
            matched.append(graph3)
        edges_combinations.clear()
        edges_combinations.extend(matched)
    
    if job_func:
        job_func("Verify the graphs...", len(edges_combinations))
    
    cdef list answer = []
    for graph1 in edges_combinations:
        if stop_func and stop_func():
            break
        if not graph1.is_connected():
            continue
        if is_isomorphic(graph1, answer):
            continue
        answer.append(graph1)
    return answer, (time() - t0)
