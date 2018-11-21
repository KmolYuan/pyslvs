# -*- coding: utf-8 -*-
# cython: language_level=3

"""Planarity check function for Graph class.

The algorithms references:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from typing import Tuple
from collections import defaultdict
from libcpp.list cimport list as c_list
from libcpp.map cimport map as c_map


cpdef inline bool is_planar(Graph g):
    """Return True if the graph is planar.
    
    Only this function well be shown as public.
    """
    cdef LRPlanarity planarity_state = LRPlanarity(g)
    return planarity_state.lr_planarity() is not None


cdef inline ConflictPair top_of_stack(list l):
    """Returns the conflict pair element on top of the stack."""
    if not l:
        return None
    return l[-1]


cdef class LRPlanarity:

    """A class to maintain the state during planarity check."""

    cdef Graph g, DG
    cdef PlanarEmbedding embedding
    cdef list S
    cdef dict lowpt, lowpt2, nesting_depth, adjs
    cdef dict ordered_adjs, stack_bottom, lowpt_edge, left_ref, right_ref

    # default dict
    cdef object height, parent_edge, ref, side

    def __cinit__(self, g: Graph):
        # copy G without adding self-loops
        self.g = g.copy()

        def return_none():
            return None

        # distance from tree root
        self.height = defaultdict(return_none)

        self.lowpt = {}  # height of lowest return point of an edge
        self.lowpt2 = {}  # height of second lowest return point
        self.nesting_depth = {}  # for nesting order

        # None -> missing edge
        self.parent_edge = defaultdict(return_none)

        # oriented DFS graph
        self.DG = Graph([])
        self.DG.add_nodes_from(g.nodes)

        self.adjs = {}
        self.ordered_adjs = {}

        self.ref = defaultdict(return_none)
        self.side = defaultdict(lambda: 1)

        # stack of conflict pairs
        self.S = []
        self.stack_bottom = {}
        self.lowpt_edge = {}

        self.left_ref = {}
        self.right_ref = {}

        self.embedding = PlanarEmbedding([])

    cdef PlanarEmbedding lr_planarity(self):
        """Execute the LR planarity test."""
        if len(self.g.nodes) > 2 and len(self.g.edges) > 3 * len(self.g.nodes) - 6:
            # graph is not planar
            return None

        # make adjacency lists for dfs
        cdef int v
        for v in self.g.nodes:
            self.adjs[v] = list(self.g.adj[v])

        # orientation of the graph by depth first search traversal
        cdef c_list[int] roots
        for v in self.g.nodes:
            if self.height[v] is None:
                self.height[v] = 0
                roots.push_back(v)
                self.dfs_orientation(v)

        # Free no longer used variables
        self.g = None
        self.lowpt2 = None
        self.adjs = None

        # testing
        cdef int n1, n2
        for v in self.DG.nodes:  # sort the adjacency lists by nesting depth
            # note: this sorting leads to non linear time
            self.ordered_adjs[v] = sorted(
                [n2 for n1, n2 in self.DG.edges if n1 == v],
                key=lambda x: self.nesting_depth[v, x]
            )
        for v in roots:
            if not self.dfs_testing(v):
                return None

        # Free no longer used variables
        self.height = None
        self.lowpt = None
        self.S = None
        self.stack_bottom = None
        self.lowpt_edge = None

        cdef tuple e
        for e in self.DG.edges:
            self.nesting_depth[e] = self.sign(e) * self.nesting_depth[e]

        self.embedding.add_nodes_from(self.DG.nodes)

        cdef int previous_node, w
        for v in self.DG.nodes:
            # sort the adjacency lists again
            self.ordered_adjs[v] = sorted(
                [n2 for n1, n2 in self.DG.edges if n1 == v],
                key=lambda x: self.nesting_depth[v, x]
            )
            # initialize the embedding
            previous_node = -1
            for w in self.ordered_adjs[v]:
                self.embedding.add_half_edge_cw(v, w, previous_node)
                previous_node = w

        # Free no longer used variables
        self.DG = None
        self.nesting_depth = None
        self.ref = None

        # compute the complete embedding
        for v in roots:
            self.dfs_embedding(v)

        return self.embedding

    cdef bool dfs_testing(self, int v):
        """Test for LR partition."""
        # the recursion stack
        cdef c_list[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef c_map[int, int] ind
        # boolean to indicate whether to skip the initial work for an edge
        cdef object skip_init = defaultdict(lambda: False)

        cdef bool skip_final
        cdef int w
        cdef tuple e, ei
        while not dfs_stack.empty():
            v = dfs_stack.back()
            dfs_stack.pop_back()
            e = self.parent_edge[v]
            # to indicate whether to skip the final block after the for loop
            skip_final = False

            for w in self.ordered_adjs[v][ind[v]:]:
                ei = (v, w)

                if not skip_init[ei]:
                    self.stack_bottom[ei] = top_of_stack(self.S)

                    if ei == self.parent_edge[w]:  # tree edge
                        dfs_stack.push_back(v)  # revisit v after finishing w
                        dfs_stack.push_back(w)  # visit w next
                        skip_init[ei] = True  # don't redo this block
                        skip_final = True  # skip final work after breaking
                        break  # handle next node in dfs_stack (i.e. w)
                    else:  # back edge
                        self.lowpt_edge[ei] = ei
                        self.S.append(ConflictPair(right=Interval(ei, ei)))

                # integrate new return edges
                if self.lowpt[ei] < self.height[v]:
                    if w == self.ordered_adjs[v][0]:  # e_i has return edge
                        self.lowpt_edge[e] = self.lowpt_edge[ei]
                    else:  # add constraints of e_i
                        if not self.add_constraints(ei, e):
                            # graph is not planar
                            return False

                ind[v] += 1

            if not skip_final:
                # remove back edges returning to parent
                if e is not None:  # v isn't root
                    self.remove_back_edges(e)

        return True

    cdef void dfs_embedding(self, int v):
        """Completes the embedding."""
        # the recursion stack
        cdef c_list[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef c_map[int, int] ind

        cdef int w
        while not dfs_stack.empty():
            v = dfs_stack.back()
            dfs_stack.pop_back()

            for w in self.ordered_adjs[v][ind[v]:]:
                ind[v] += 1
                ei = (v, w)

                if ei == self.parent_edge[w]:  # tree edge
                    self.embedding.add_half_edge_first(w, v)
                    self.left_ref[v] = w
                    self.right_ref[v] = w

                    dfs_stack.push_back(v)  # revisit v after finishing w
                    dfs_stack.push_back(w)  # visit w next
                    break  # handle next node in dfs_stack (i.e. w)
                else:  # back edge
                    if self.side[ei] == 1:
                        self.embedding.add_half_edge_cw(w, v, self.right_ref[w])
                    else:
                        self.embedding.add_half_edge_ccw(w, v, self.left_ref[w])
                        self.left_ref[w] = v

    cdef void dfs_orientation(self, int v):
        """Orient the graph by DFS, compute lowpoints and nesting order."""
        # the recursion stack
        cdef c_list[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef c_map[int, int] ind
        # boolean to indicate whether to skip the initial work for an edge
        cdef object skip_init = defaultdict(lambda: False)

        cdef int w
        cdef tuple e, vw
        while not dfs_stack.empty():
            v = dfs_stack.back()
            dfs_stack.pop_back()
            e = self.parent_edge[v]

            for w in self.adjs[v][ind[v]:]:
                vw = (v, w)

                if not skip_init[vw]:
                    if (v, w) in self.DG.edges or (w, v) in self.DG.edges:
                        ind[v] += 1
                        continue  # the edge was already oriented

                    self.DG.add_edge(v, w)  # orient the edge

                    self.lowpt[vw] = self.height[v]
                    self.lowpt2[vw] = self.height[v]
                    if self.height[w] is None:  # (v, w) is a tree edge
                        self.parent_edge[w] = vw
                        self.height[w] = self.height[v] + 1

                        dfs_stack.push_back(v)  # revisit v after finishing w
                        dfs_stack.push_back(w)  # visit w next
                        skip_init[vw] = True  # don't redo this block
                        break  # handle next node in dfs_stack (i.e. w)
                    else:  # (v, w) is a back edge
                        self.lowpt[vw] = self.height[w]

                # determine nesting graph
                self.nesting_depth[vw] = 2 * self.lowpt[vw]
                if self.lowpt2[vw] < self.height[v]:  # chordal
                    self.nesting_depth[vw] += 1

                # update lowpoints of parent edge e
                if e is not None:
                    if self.lowpt[vw] < self.lowpt[e]:
                        self.lowpt2[e] = min(self.lowpt[e], self.lowpt2[vw])
                        self.lowpt[e] = self.lowpt[vw]
                    elif self.lowpt[vw] > self.lowpt[e]:
                        self.lowpt2[e] = min(self.lowpt2[e], self.lowpt[vw])
                    else:
                        self.lowpt2[e] = min(self.lowpt2[e], self.lowpt2[vw])

                ind[v] += 1

    cdef bool add_constraints(self, tuple ei, tuple e):
        cdef ConflictPair P = ConflictPair()

        # merge return edges of e_i into P.right
        cdef ConflictPair Q
        while True:
            Q = self.S.pop()
            if not Q.left.empty():
                Q.swap()
            if not Q.left.empty():  # not planar
                return False
            if self.lowpt[Q.right.low] > self.lowpt[e]:
                # merge intervals
                if P.right.empty():  # topmost interval
                    P.right = Q.right.copy()
                else:
                    self.ref[P.right.low] = Q.right.high
                P.right.low = Q.right.low
            else:  # align
                self.ref[Q.right.low] = self.lowpt_edge[e]
            if top_of_stack(self.S) == self.stack_bottom[ei]:
                break

        # merge conflicting return edges of e_1,...,e_i-1 into P.L
        while (
            top_of_stack(self.S).left.conflicting(ei, self) or
            top_of_stack(self.S).right.conflicting(ei, self)
        ):
            Q = self.S.pop()
            if Q.right.conflicting(ei, self):
                Q.swap()
            if Q.right.conflicting(ei, self):  # not planar
                return False
            # merge interval below lowpt(e_i) into P.R
            self.ref[P.right.low] = Q.right.high
            if Q.right.low is not None:
                P.right.low = Q.right.low

            if P.left.empty():  # topmost interval
                P.left = Q.left.copy()
            else:
                self.ref[P.left.low] = Q.left.high
            P.left.low = Q.left.low

        if not (P.left.empty() and P.right.empty()):
            self.S.append(P)

        return True

    cdef void remove_back_edges(self, tuple e):
        cdef int u = e[0]

        # trim back edges ending at parent u
        # drop entire conflict pairs
        cdef ConflictPair p
        while self.S and top_of_stack(self.S).lowest(self) == self.height[u]:
            p = self.S.pop()
            if p.left.low is not None:
                self.side[p.left.low] = -1

        if self.S:  # one more conflict pair to consider
            p = self.S.pop()
            # trim left interval
            while p.left.high is not None and p.left.high[1] == u:
                p.left.high = self.ref[p.left.high]
            if p.left.high is None and p.left.low is not None:
                # just emptied
                self.ref[p.left.low] = p.right.low
                self.side[p.left.low] = -1
                p.left.low = None
            # trim right interval
            while p.right.high is not None and p.right.high[1] == u:
                p.right.high = self.ref[p.right.high]
            if p.right.high is None and p.right.low is not None:
                # just emptied
                self.ref[p.right.low] = p.left.low
                self.side[p.right.low] = -1
                p.right.low = None
            self.S.append(p)

        # side of e is side of a highest return edge
        cdef tuple hl, hr
        if self.lowpt[e] < self.height[u]:  # e has return edge
            p = top_of_stack(self.S)
            hl = p.left.high
            hr = p.right.high

            if (hl is not None) and (hr is None or self.lowpt[hl] > self.lowpt[hr]):
                self.ref[e] = hl
            else:
                self.ref[e] = hr

    cdef int sign(self, tuple e):
        """Resolve the relative side of an edge to the absolute side."""
        # the recursion stack
        cdef list dfs_stack = [e]
        # dict to remember reference edges
        cdef object old_ref = defaultdict(lambda: None)

        while dfs_stack:
            e = dfs_stack.pop()

            if self.ref[e] is not None:
                dfs_stack.append(e)  # revisit e after finishing self.ref[e]
                dfs_stack.append(self.ref[e])  # visit self.ref[e] next
                old_ref[e] = self.ref[e]  # remember value of self.ref[e]
                self.ref[e] = None
            else:
                self.side[e] *= self.side[old_ref[e]]

        return self.side[e]


cdef class ConflictPair:

    """Represents a different constraint between two intervals.

    The edges in the left interval must have a different orientation than
    the one in the right interval.
    """

    cdef Interval left, right

    def __cinit__(self, left=Interval(), right=Interval()):
        self.left = left
        self.right = right

    cdef void swap(self):
        """Swap left and right intervals"""
        self.left, self.right = self.right, self.left

    cdef int lowest(self, LRPlanarity planarity_state):
        """Return the lowest low point of a conflict pair"""
        if self.left.empty():
            return planarity_state.lowpt[self.right.low]
        if self.right.empty():
            return planarity_state.lowpt[self.left.low]
        return min(
            planarity_state.lowpt[self.left.low],
            planarity_state.lowpt[self.right.low]
        )


cdef class Interval:

    """Represents a set of return edges.

    All return edges in an interval induce a same constraint on the contained
    edges, which means that all edges must either have a left orientation or
    all edges must have a right orientation.
    """

    cdef tuple low, high

    def __cinit__(self, low: Tuple[int, int] = None, high: Tuple[int, int] = None):
        self.low = low
        self.high = high

    cdef bool empty(self):
        """Check if the interval is empty"""
        return (self.low is None) and (self.high is None)

    cdef Interval copy(self):
        """Return a copy of this interval"""
        return Interval(self.low, self.high)

    cdef bool conflicting(self, tuple b, LRPlanarity state):
        """Return True if interval I conflicts with edge b"""
        return not self.empty() and state.lowpt[self.high] > state.lowpt[b]


cdef class PlanarEmbedding(Graph):

    """Represents a planar graph with its planar embedding."""

    cdef c_map[int, int] node_label
    cdef c_map[int, c_map[int, c_map[int, int]]] edge_label

    cdef void add_half_edge_cw(self, int start_node, int end_node, int reference_neighbor):
        """Adds a half-edge from start_node to end_node.
        
        edge_label[:, :, :] 0: cw
        edge_label[:, :, :] 1: ccw
        node_label[:]: first_nbr
        """
        self.add_edge(start_node, end_node)  # Add edge to graph

        if reference_neighbor == -1:
            # The start node has no neighbors
            self.edge_label[start_node][end_node][0] = end_node
            self.edge_label[start_node][end_node][1] = end_node
            self.node_label[start_node] = end_node
            return

        # Get half-edge at the other side
        cdef int cw_reference = self.edge_label[start_node][reference_neighbor][0]
        # Alter half-edge data structures
        self.edge_label[start_node][reference_neighbor][0] = end_node
        self.edge_label[start_node][end_node][0] = cw_reference
        self.edge_label[start_node][cw_reference][1] = end_node
        self.edge_label[start_node][end_node][1] = reference_neighbor

    cdef void add_half_edge_ccw(self, int start_node, int end_node, int reference_neighbor):
        """Adds a half-edge from start_node to end_node."""
        self.add_half_edge_cw(
            start_node,
            end_node,
            self.edge_label[start_node][reference_neighbor][1]
        )
        cdef int s_f_label
        if self.node_label.find(start_node) != self.node_label.end():
            s_f_label = self.node_label[start_node]
        else:
            s_f_label = -1
        if reference_neighbor == s_f_label:
            # Update first neighbor
            self.node_label[start_node] = end_node

    cdef void add_half_edge_first(self, int start_node, int end_node):
        """The added half-edge is inserted at the first position in the order."""
        self.add_half_edge_ccw(start_node, end_node, self.node_label[start_node])
