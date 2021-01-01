# -*- coding: utf-8 -*-
# cython: language_level=3

"""Planarity check function for Graph class.

The algorithms references:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from collections import defaultdict
from libcpp.pair cimport pair
from libcpp.list cimport list as clist
from libcpp.map cimport map

ctypedef pair[int, int] ipair


cpdef bint is_planar(Graph g):
    """Return true if the graph is a planar graph."""
    cdef _LRPlanarity planarity_state = _LRPlanarity.__new__(_LRPlanarity, g)
    return planarity_state.lr_planarity() is not None


cdef inline _ConflictPair _stack_top(list l):
    """Returns the conflict pair element on top of the stack."""
    if not l:
        return None
    return l[-1]


cdef class _LRPlanarity:
    """A class to maintain the state during planarity check."""
    cdef Graph g, DG
    cdef _PlanarEmbedding embedding
    cdef list s
    cdef dict lowpt, lowpt2, nesting_depth, adjs
    cdef dict ordered_adjs, stack_bottom, lowpt_edge, left_ref, right_ref
    # default dict
    cdef object height, parent_edge, ref, side

    def __cinit__(self, Graph g):
        # copy G without adding self-loops
        self.g = Graph.__new__(Graph, set(g.edges))
        # distance from tree root
        self.height = defaultdict(lambda: None)
        self.lowpt = {}  # height of lowest return point of an edge
        self.lowpt2 = {}  # height of second lowest return point
        self.nesting_depth = {}  # for nesting order
        # None -> missing edge
        self.parent_edge = defaultdict(lambda: None)
        # oriented DFS graph
        self.DG = Graph.__new__(Graph, [])
        self.DG.add_vertices(g.vertices)
        self.adjs = {}
        self.ordered_adjs = {}
        self.ref = defaultdict(lambda: None)
        self.side = defaultdict(lambda: 1)
        # stack of conflict pairs
        self.s = []
        self.stack_bottom = {}
        self.lowpt_edge = {}
        self.left_ref = {}
        self.right_ref = {}
        self.embedding = _PlanarEmbedding.__new__(_PlanarEmbedding, [])

    cdef _PlanarEmbedding lr_planarity(self):
        """Execute the LR planarity test."""
        if len(self.g.vertices) > 2 and len(self.g.edges) > 3 * len(self.g.vertices) - 6:
            # graph is not planar
            return None
        # make adjacency lists for dfs
        cdef int v
        for v in self.g.vertices:
            self.adjs[v] = list(self.g.adj[v])
        # orientation of the graph by depth first search traversal
        cdef clist[int] roots
        for v in self.g.vertices:
            if self.height[v] is None:
                self.height[v] = 0
                roots.push_back(v)
                self.dfs_orientation(v)
        # Free no longer used variables
        self.g = None
        self.lowpt2 = None
        self.adjs = None
        # testing
        for v in self.DG.vertices:  # sort the adjacency lists by nesting depth
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
        self.s = None
        self.stack_bottom = None
        self.lowpt_edge = None
        for e in self.DG.edges:
            self.nesting_depth[e] = self.sign(e) * self.nesting_depth[e]
        self.embedding.add_vertices(self.DG.vertices)
        cdef int previous_node, w
        for v in self.DG.vertices:
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

    cdef bint dfs_testing(self, int v):
        """Test for LR partition."""
        # the recursion stack
        cdef clist[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef map[int, int] ind
        # bintean to indicate whether to skip the initial work for an edge
        skip_init = defaultdict(lambda: False)
        cdef bint skip_final
        cdef int w
        while not dfs_stack.empty():
            v = dfs_stack.back()
            dfs_stack.pop_back()
            e = self.parent_edge[v]
            # to indicate whether to skip the final block after the for loop
            skip_final = False
            for w in self.ordered_adjs[v][ind[v]:]:
                ei = (v, w)
                if not skip_init[ei]:
                    self.stack_bottom[ei] = _stack_top(self.s)
                    if ei == self.parent_edge[w]:  # tree edge
                        dfs_stack.push_back(v)  # revisit v after finishing w
                        dfs_stack.push_back(w)  # visit w next
                        skip_init[ei] = True  # don't redo this block
                        skip_final = True  # skip final work after breaking
                        break  # handle next node in dfs_stack (i.e. w)
                    else:  # back edge
                        self.lowpt_edge[ei] = ei
                        self.s.append(_ConflictPair.__new__(_ConflictPair, None,
                            _Interval.__new__(_Interval, ei, ei)))
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
        cdef clist[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef map[int, int] ind
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
        cdef clist[int] dfs_stack = [v]
        # index of next edge to handle in adjacency list of each node
        cdef map[int, int] ind
        # bintean to indicate whether to skip the initial work for an edge
        skip_init = defaultdict(lambda: False)
        cdef int w
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

    cdef bint add_constraints(self, tuple ei, tuple e):
        cdef _ConflictPair p = _ConflictPair.__new__(_ConflictPair, None, None)
        # merge return edges of e_i into P.right
        cdef _ConflictPair q
        while True:
            q = self.s.pop()
            if not q.left.empty():
                q.swap()
            if not q.left.empty():  # not planar
                return False
            if self.lowpt[q.right.low] > self.lowpt[e]:
                # merge intervals
                if p.right.empty():  # topmost interval
                    p.right = q.right.copy()
                else:
                    self.ref[p.right.low] = q.right.high
                p.right.low = q.right.low
            else:  # align
                self.ref[q.right.low] = self.lowpt_edge[e]
            if _stack_top(self.s) == self.stack_bottom[ei]:
                break
        # merge conflicting return edges of e_1,...,e_i-1 into P.L
        while (
            _stack_top(self.s).left.conflicting(ei, self) or
            _stack_top(self.s).right.conflicting(ei, self)
        ):
            q = self.s.pop()
            if q.right.conflicting(ei, self):
                q.swap()
            if q.right.conflicting(ei, self):  # not planar
                return False
            # merge interval below lowpt(e_i) into P.R
            self.ref[p.right.low] = q.right.high
            if q.right.low is not None:
                p.right.low = q.right.low
            if p.left.empty():  # topmost interval
                p.left = q.left.copy()
            else:
                self.ref[p.left.low] = q.left.high
            p.left.low = q.left.low
        if not (p.left.empty() and p.right.empty()):
            self.s.append(p)
        return True

    cdef void _trim_interval(self, _Interval inter1, _Interval inter2, int u):
        """Trim left or right interval."""
        while inter1.high is not None and inter1.high[1] == u:
            inter1.high = self.ref[inter1.high]
        if inter1.high is None and inter1.low is not None:
            # just emptied
            self.ref[inter1.low] = inter2.low
            self.side[inter1.low] = -1
            inter1.low = None

    cdef void remove_back_edges(self, tuple e):
        cdef int u = e[0]
        # trim back edges ending at parent u
        # drop entire conflict pairs
        cdef _ConflictPair p
        while self.s and _stack_top(self.s).lowest(self) == self.height[u]:
            p = self.s.pop()
            if p.left.low is not None:
                self.side[p.left.low] = -1
        if self.s:  # one more conflict pair to consider
            p = self.s.pop()
            # trim left interval
            self._trim_interval(p.left, p.right, u)
            # trim right interval
            self._trim_interval(p.right, p.left, u)
            self.s.append(p)
        # side of e is side of a highest return edge
        if self.lowpt[e] < self.height[u]:  # e has return edge
            p = _stack_top(self.s)
            hl = p.left.high
            hr = p.right.high
            if (hl is not None) and (hr is None or self.lowpt[hl] > self.lowpt[hr]):
                self.ref[e] = hl
            else:
                self.ref[e] = hr

    cdef int sign(self, tuple e):
        """Resolve the relative side of an edge to the absolute side."""
        # the recursion stack
        dfs_stack = [e]
        # dict to remember reference edges
        old_ref = defaultdict(lambda: None)
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


cdef class _ConflictPair:
    """Represents a different constraint between two intervals.

    The edges in the left interval must have a different orientation than
    the one in the right interval.
    """
    cdef _Interval left, right

    def __cinit__(self, _Interval left, _Interval right):
        if left is None:
            self.left = _Interval.__new__(_Interval, None, None)
        else:
            self.left = left
        if right is None:
            self.right = _Interval.__new__(_Interval, None, None)
        else:
            self.right = right

    cdef void swap(self):
        """Swap left and right intervals"""
        self.left, self.right = self.right, self.left

    cdef int lowest(self, _LRPlanarity planarity_state):
        """Return the lowest low point of a conflict pair"""
        if self.left.empty():
            return planarity_state.lowpt[self.right.low]
        if self.right.empty():
            return planarity_state.lowpt[self.left.low]
        return min(
            planarity_state.lowpt[self.left.low],
            planarity_state.lowpt[self.right.low]
        )


cdef class _Interval:
    """Represents a set of return edges.

    All return edges in an interval induce a same constraint on the contained
    edges, which means that all edges must either have a left orientation or
    all edges must have a right orientation.
    """
    cdef tuple low, high

    def __cinit__(self, tuple low, tuple high):
        self.low = low
        self.high = high

    cdef bint empty(self):
        """Check if the interval is empty"""
        return self.low is self.high is None

    cdef _Interval copy(self):
        """Return a copy of this interval"""
        return _Interval.__new__(_Interval, self.low, self.high)

    cdef bint conflicting(self, tuple b, _LRPlanarity state):
        """Return true if interval I conflicts with edge b"""
        return not self.empty() and state.lowpt[self.high] > state.lowpt[b]


cdef class _PlanarEmbedding(Graph):
    """Represents a planar graph with its planar embedding."""
    cdef map[int, int] node_label
    cdef map[ipair, ipair] edge_label

    cdef void add_half_edge_cw(self, int start_node, int end_node, int reference_neighbor):
        """Adds a half-edge from start_node to end_node.

        edge_label[:, :, :] 0: cw
        edge_label[:, :, :] 1: ccw
        node_label[:]: first_nbr
        """
        self.add_edge(start_node, end_node)  # Add edge to graph
        if reference_neighbor == -1:
            # The start node has no neighbors
            self.edge_label[ipair(start_node, end_node)].first = end_node
            self.edge_label[ipair(start_node, end_node)].second = end_node
            self.node_label[start_node] = end_node
            return
        # Get half-edge at the other side
        cdef int cw_reference = self.edge_label[ipair(start_node, reference_neighbor)].first
        # Alter half-edge data structures
        self.edge_label[ipair(start_node, reference_neighbor)].first = end_node
        self.edge_label[ipair(start_node, end_node)].first = cw_reference
        self.edge_label[ipair(start_node, cw_reference)].second = end_node
        self.edge_label[ipair(start_node, end_node)].second = reference_neighbor

    cdef void add_half_edge_ccw(self, int start_node, int end_node, int reference_neighbor):
        """Adds a half-edge from start_node to end_node."""
        self.add_half_edge_cw(
            start_node,
            end_node,
            self.edge_label[ipair(start_node, reference_neighbor)].second
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
