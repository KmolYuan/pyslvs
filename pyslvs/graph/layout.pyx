# -*- coding: utf-8 -*-
# cython: language_level=3

"""Layout functions for graph object.

The algorithms references:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from typing import Set, Iterable, Any
cimport cython
from cpython cimport PyDict_Contains, PyIndex_Check
from cpython.slice cimport PySlice_GetIndicesEx
from libc.math cimport hypot, M_PI, sin, cos, atan2
from libcpp.map cimport map
from .graph cimport Graph


cpdef dict external_loop_layout(Graph g, bint node_mode, double scale = 1.):
    """Layout position decided by outer loop (max cycle).

    Return the layout position decided by external loop.
    Argument `node_mode` will transform edges into vertices.
    Argument `scale` will resize the position by scale factor.
    """
    pos = _line_polygon_layout(g, scale)
    # Node mode
    cdef int i, n1, n2
    cdef double x1, y1, x2, y2
    if not node_mode:
        n_pos = {}
        edges_view = sorted([((n2, n1) if n2 < n1 else (n1, n2)) for (n1, n2) in g.edges])
        for i, (n1, n2) in enumerate(edges_view):
            x1, y1 = pos[n1]
            x2, y2 = pos[n2]
            n_pos[i] = _middle_point(x1, y1, x2, y2)
        pos = n_pos

    return pos


cdef inline dict _line_polygon_layout(Graph g, double scale):
    """Generate position of external loop and inner lines."""
    cdef OrderedSet o_loop = _external_loop(g)
    lines = None
    while lines is None:
        # Patch function of external cycle
        lines = _inner_lines(g, o_loop)

    o_loop.roll(min(o_loop), 0)
    pos = {}
    _regular_polygon_layout(o_loop, scale, pos)
    cdef OrderedSet used_nodes = o_loop.copy()

    cdef int start, end
    cdef OrderedSet line
    cdef map[int, map[int, int]] line_limit, line_counter
    for line, start, end in lines:
        line_limit[start][end] += 1

    cdef int limit
    cdef double p, x1, y1, x2, y2
    for line, start, end in lines:
        x1, y1 = pos[start]
        x2, y2 = pos[end]
        # Conditions of linear or bezier curve
        limit = line_limit[start][end]
        if limit == 1:
            _linear_layout(x1, y1, x2, y2, line, pos)
        else:
            line_counter[start][end] += 1
            p = <double>line_counter[start][end] / (limit + 1)
            _bezier_layout(p, x1, y1, x2, y2, line, pos)
        used_nodes.update(line)

    # Last check for debug
    if set(g.vertices) != set(pos):
        raise ValueError(
            f"the algorithm is error with {g.edges}\n"
            f"external loop: {o_loop}\n"
            f"inner layout: {set(pos) - o_loop}\n"
            f"node {set(g.vertices) - set(pos)} are not included"
        )

    return pos


cdef inline void _regular_polygon_layout(OrderedSet vertices, double scale, dict pos):
    """Return position of a regular polygon with radius 5.
    Start from bottom with clockwise.
    """
    cdef int edge_count = len(vertices)
    scale *= 5
    cdef int i
    cdef double angle = M_PI * 1.5
    cdef double angle_step = M_PI * 2 / edge_count
    for i in range(edge_count):
        pos[vertices[i]] = (scale * cos(angle), scale * sin(angle))
        angle -= angle_step


cdef inline void _linear_layout(
    double x1,
    double y1,
    double x2,
    double y2,
    OrderedSet vertices,
    dict pos
):
    """Layout position decided by equal division between two points."""
    cdef int count = len(vertices)
    if count < 1:
        raise ValueError(f"Invalid point number {count}")

    count += 1
    cdef double sx = (x2 - x1) / count
    cdef double sy = (y2 - y1) / count

    cdef int i
    for i in range(1, count):
        pos[vertices[i - 1]] = (x1 + i * sx, y1 + i * sy)


cdef inline void _bezier_layout(
    double p,
    double x1,
    double y1,
    double x2,
    double y2,
    OrderedSet vertices,
    dict pos
):
    """Layout position decided by four points bezier curve."""
    cdef int count = len(vertices)
    if count < 1:
        raise ValueError(f"Invalid point number {count}")
    count += 1

    # Make direction from start and end nodes
    cdef double sx = x2 - x1
    cdef double sy = y2 - y1
    cdef double r = hypot(sx, sy) * 0.45
    cdef double direction = (p + 0.5) * M_PI + atan2(sy, sx)
    cdef double cx1 = x1 + r * cos(direction)
    cdef double cy1 = y1 + r * sin(direction)
    cdef double cx2 = x2 + r * cos(-direction)
    cdef double cy2 = y2 + r * sin(-direction)

    cdef int i
    cdef double u
    for i in range(1, count):
        u = <double>i / count
        pos[vertices[i - 1]] = _bezier_curve(u, x1, y1, cx1, cy1, cx2, cy2, x2, y2)


cdef inline tuple _bezier_curve(
    double u,
    double x1,
    double y1,
    double x2,
    double y2,
    double x3,
    double y3,
    double x4,
    double y4
):
    """A coordinate on bezier curve."""
    cdef double mx = _uv(u, x2, x3)
    cdef double my = _uv(u, y2, y3)
    return (
        _uv(u, _uv(u, _uv(u, x1, x2), mx), _uv(u, mx, _uv(u, x3, x4))),
        _uv(u, _uv(u, _uv(u, y1, y2), my), _uv(u, my, _uv(u, y3, y4)))
    )


cdef inline double _uv(double u, double v1, double v2) nogil:
    """Return v1 + u * (v2 - v1)"""
    return v1 + u * (v2 - v1)


cdef inline tuple _middle_point(double x1, double y1, double x2, double y2):
    """Return middle point of two coordinates."""
    return ((x1 + x2) / 2), ((y1 + y2) / 2)


cdef list _inner_lines(Graph g, OrderedSet o_loop):
    """Layout for inner vertices of graph block."""
    cdef OrderedSet vertices = OrderedSet(g.vertices) - o_loop
    if not vertices:
        return []
    lines = []
    cdef OrderedSet used_nodes = o_loop.copy()
    cdef int n
    cdef OrderedSet line, inter
    while vertices:
        n = vertices.pop(0)
        if not (used_nodes & g.adj[n]):
            # Not contacted yet
            vertices.add(n)
            continue

        line = OrderedSet.__new__(OrderedSet)
        line.add(n)
        inter = OrderedSet(g.adj[n]) - used_nodes
        while inter:
            # New vertices to add
            n = inter.pop()
            line.add(n)
            vertices.remove(n)
            inter = vertices & g.adj[n]
            if used_nodes & g.adj[n]:
                # Connected with any used node
                break

        # Find the intersections of the line
        if line[0] == line[-1]:
            inter = used_nodes & g.adj[line[0]]
            start = inter.pop()
            end = inter.pop()
        else:
            inter = used_nodes & g.adj[line[0]]
            start = inter.pop()
            inter = used_nodes & g.adj[line[-1]]
            end = inter.pop()

        if _split_loop(o_loop, line, start, end):
            # 'o_loop' has been changed
            return None

        # A unified format of start node and end node
        if end < start:
            end, start = start, end
            line.reverse()

        # Line is ended
        used_nodes.update(line)
        lines.append((line, start, end))

    return lines


cdef inline bint _split_loop(OrderedSet o_loop, OrderedSet line, int n1, int n2):
    """Return false if the line is a chord.
    Or it is an arc of external cycle.
    """
    if n1 == n2:
        return False
    loop_list = list(o_loop) * 2
    cdef int i, n
    cdef int s0 = -1
    cdef int s1 = -1
    cdef int s2 = -1
    for i, n in enumerate(loop_list):
        if n not in {n1, n2}:
            continue
        if s0 == -1:
            s0 = i
        elif s1 == -1:
            s1 = i
        else:
            s2 = i
            break
    else:
        # Split failed
        return False

    # Remove the last repeated part
    loop_list = loop_list[:s2 + 1]

    cdef int anchor
    if len(line) > s1 - (s0 + 1):
        anchor = s0 + 1
        del loop_list[anchor:s1]
        if loop_list[s0] == n2:
            # reverse
            loop_list[anchor:anchor] = reversed(line)
        else:
            loop_list[anchor:anchor] = line
    elif len(line) > s2 - (s1 + 1):
        anchor = s1 + 1
        del loop_list[anchor:s2]
        if loop_list[s1] == n2:
            # reverse
            loop_list[anchor:anchor] = reversed(line)
        else:
            loop_list[anchor:anchor] = line
    else:
        # Is a normal chord
        return False

    o_loop.clear()
    # Remove the first repeated part
    o_loop.update(loop_list[s0:])
    return True


cdef inline OrderedSet _external_loop(Graph g):
    """Return vertices of external loop."""
    cycles = _cycle_basis(g)
    if not cycles:
        raise ValueError(f"invalid graph has no any cycle: {g.edges}")

    cdef OrderedSet c1, c2
    while len(cycles) > 1:
        c1 = cycles.pop()
        for c2 in cycles:
            if _merge_inter(c1, c2, cycles, g):
                break
        else:
            # Cycles have no contacted
            # Find connection from edges
            for c2 in cycles:
                if _merge_no_inter(c1, c2, cycles, g):
                    break
            else:
                raise ValueError(
                    f"invalid graph: {g.edges}\n"
                    f"last one: {c1}\n"
                    f"with cycle(s): {cycles}"
                )
    return cycles.pop()


cdef inline bint _merge_inter(
    OrderedSet c1,
    OrderedSet c2,
    list cycles,
    Graph g
):
    """Merge function for intersection strategy."""
    # Find the intersection
    if not len(c1 & c2) >= 2:
        return False

    # Ignore subsets
    if c1 >= c2:
        cycles.remove(c2)
        cycles.append(c1)
        return True
    inter_over = c1.ordered_intersections(c2, is_loop=True)
    # Find the intersection with reversed cycle
    c2.reverse()
    inter_tmp = c1.ordered_intersections(c2, is_loop=True)
    if len(inter_tmp) < len(inter_over):
        # Choose the longest continuous intersection
        inter_over = inter_tmp
    inter_tmp = None

    cdef int start = -1
    cdef int end = -1

    cdef int i
    cdef OrderedSet inter
    for i, inter in enumerate(inter_over):
        if not inter.is_ordered_subset(c1, is_loop=True):
            # Intersection and cycle 1 has wrong direction
            inter.reverse()
        if inter.is_ordered_subset(c2, is_loop=True) and i == 0:
            # Cycle 1 and cycle 2 should has different direction
            c2.reverse()

        # Prune cycle 2 by intersection
        c2 -= inter[1:-1]

        # Interface nodes
        if i == 0:
            start = inter[0]
        end = inter[-1]

    # Roll to interface
    c1.roll(end, -1)
    c2.roll(start, 0)

    # Insert points
    _compare_insert(
        c1,
        c2,
        c1.index(start),
        c1.index(end),
        c2.index(start),
        c2.index(end),
        g
    )

    # The cycle 2 has been merged into cycle 1
    cycles.remove(c2)
    cycles.append(c1)
    return True


cdef inline bint _merge_no_inter(
    OrderedSet c1,
    OrderedSet c2,
    list cycles,
    Graph g
):
    """Merge function for the strategy without intersection."""
    cdef OrderedSet inter = OrderedSet.__new__(OrderedSet)
    cdef map[int, int] inter_map

    cdef int n1, n2
    for n1, n2 in g.edges:
        if (n1 in c1) and (n2 in c2):
            inter.add(n1)
            inter_map[n1] = n2
        elif (n2 in c1) and (n1 in c2):
            inter.add(n2)
            inter_map[n2] = n1

    if not inter:
        return False

    # Resort intersection
    if not inter.is_ordered_subset(c1, is_loop=True):
        inter = c1 & inter
    cdef int start = inter[0]
    cdef int end = inter[-1]
    cdef int replace_start = c2.index(inter_map[start])
    cdef int replace_end = c2.index(inter_map[end])
    if replace_start > replace_end:
        c2.reverse()

    # Roll to interface
    c1.roll(end, -1)
    c2.roll(inter_map[start], 0)

    # Merge them
    _compare_insert(
        c1,
        c2,
        c1.index(start),
        c1.index(end),
        replace_start,
        replace_end,
        g
    )

    # The cycle 2 has been merged into cycle 1
    cycles.remove(c2)
    cycles.append(c1)
    return True


cdef inline void _compare_insert(
    OrderedSet c1,
    OrderedSet c2,
    int insert_start,
    int insert_end,
    int replace_start,
    int replace_end,
    Graph g
):
    """Compare and insert cycle 2 to cycle 1."""
    cdef int c1_degrees = 0
    cdef int c2_degrees = 0
    cdef OrderedSet c1_slice = c1[insert_start:insert_end]
    cdef OrderedSet c2_slice = c2[replace_start:replace_end]

    cdef int i
    other_nodes = set(c1_slice | c2_slice)
    for i in c1_slice:
        c1_degrees += len(set(g.adj[i]) - other_nodes)
    for i in c2_slice:
        c2_degrees += len(set(g.adj[i]) - other_nodes)

    if not (c2_degrees > c1_degrees):
        return

    # Cycle 2 should longer then intersection
    del c1[insert_start:insert_end]
    c1.insert_from(insert_start, c2_slice)


cdef inline list _cycle_basis(Graph g):
    """Returns a list of cycles which form a basis for cycles of G.
    Reference from NetworkX.
    """
    g_nodes = set(g.vertices)
    cycles = []
    cdef int root = -1
    cdef int z, nbr, p
    cdef OrderedSet cycle
    while g_nodes:
        # loop over connected components
        if root == -1:
            root = g_nodes.pop()
        stack = [root]
        pred = {root: root}
        used = {root: set()}
        # walk the spanning tree finding cycles
        while stack:
            # use last-in so cycles easier to find
            z = stack.pop()
            zused = used[z]
            for nbr in g.adj[z]:
                if nbr not in used:
                    # new node
                    pred[nbr] = z
                    stack.append(nbr)
                    used[nbr] = {z}
                elif nbr == z:
                    # self loops
                    cycles.append(OrderedSet([z]))
                elif nbr not in zused:
                    # found a cycle
                    pn = used[nbr]
                    cycle = OrderedSet([nbr, z])
                    p = pred[z]
                    while p not in pn:
                        cycle.add(p)
                        p = pred[p]
                    cycle.add(p)
                    cycles.append(cycle)
                    used[nbr].add(z)
        g_nodes -= set(pred)
        root = -1
    return cycles


cdef class _Entry:
    cdef object key
    cdef _Entry prev
    cdef _Entry next


cdef inline void _add(OrderedSet oset, object key):
    if PyDict_Contains(oset.map, key):
        return
    cdef _Entry next_entry = _Entry()
    next_entry.key = key
    next_entry.prev = oset.end.prev
    next_entry.next = oset.end
    oset.end.prev.next = oset.end.prev = oset.map[key] = next_entry
    oset.os_used += 1


cdef void _discard(OrderedSet oset, object key):
    if not PyDict_Contains(oset.map, key):
        return
    cdef _Entry entry = oset.map.pop(key)
    entry.prev.next = entry.next
    entry.next.prev = entry.prev
    oset.os_used -= 1


cdef inline bint _isorderedsubset(seq1, seq2, bint is_loop):
    """Return true if 'seq1' is ordered subset of 'seq2'."""
    cdef int seq1_len = len(seq1)
    if not seq1_len <= len(seq2):
        # 'seq1' is obviously not a subset.
        return False
    seq1_list = list(seq1)
    seq2_list = list(seq2)
    if is_loop:
        seq2_list *= 2
    cdef int matched = 0
    for self_elem in seq2_list:
        if self_elem == seq1_list[matched]:
            matched += 1
            if matched == seq1_len:
                return True
        else:
            matched = 0

    return matched == seq1_len


cdef inline bint _not_subset(OrderedSet o_set, list members):
    """Return true if 'o_set' is not any subset of members."""
    cdef OrderedSet member
    for member in members:
        if member < o_set:
            # Is superset
            members.remove(member)
            return True
        if member >= o_set:
            return False
    # No subset
    return True


cdef class _OrderedSetIterator:
    """Ordered set iterator."""
    cdef OrderedSet oset
    cdef _Entry curr
    cdef ssize_t si_used

    def __cinit__(self, OrderedSet oset):
        self.oset = oset
        self.curr = oset.end
        self.si_used = oset.os_used

    def __iter__(self):
        return self

    def __next__(self):
        cdef _Entry item

        if self.si_used != self.oset.os_used:
            # make this state sticky
            self.si_used = -1
            raise RuntimeError(f'{type(self.oset).__name__} changed size during iteration')

        item = self.curr.next
        if item == self.oset.end:
            raise StopIteration
        self.curr = item
        return item.key


cdef class _OrderedSetReverseIterator:
    """Ordered set iterator with reversed order."""
    cdef OrderedSet oset
    cdef _Entry curr
    cdef ssize_t si_used

    def __cinit__(self, OrderedSet oset):
        self.oset = oset
        self.curr = oset.end
        self.si_used = oset.os_used

    def __iter__(self):
        return self

    def __next__(self):
        if self.si_used != self.oset.os_used:
            # make this state sticky
            self.si_used = -1
            raise RuntimeError(f'{type(self.oset).__name__} changed size during iteration')

        cdef _Entry item = self.curr.prev
        if item is self.oset.end:
            raise StopIteration
        self.curr = item
        return item.key


@cython.final
cdef class OrderedSet:
    """Ordered set container."""
    cdef dict map
    cdef _Entry end
    cdef ssize_t os_used

    def __cinit__(self):
        self.map = {}
        self.os_used = 0
        self.end = end = _Entry()
        end.prev = end.next = end

    def __init__(self, iterable: Iterable[Any] = None):
        if iterable is None:
            return
        map_d = self.map
        cdef _Entry next_e
        cdef _Entry end = self.end
        for elem in iterable:
            if not PyDict_Contains(map_d, elem):
                next_e = _Entry()
                next_e.key, next_e.prev, next_e.next = elem, end.prev, end
                end.prev.next = end.prev = map_d[elem] = next_e
                self.os_used += 1

    @classmethod
    def _from_iterable(cls, it: Iterable[Any]) -> OrderedSet:
        if isinstance(it, OrderedSet):
            return it
        else:
            return cls(it)

    ##
    # set methods
    ##
    cpdef void add(self, elem):
        """Add 'elem' to the set."""
        _add(self, elem)

    cpdef void insert(self, int index, elem):
        """Insert 'elem' to 'index'."""
        self.insert_from(index, (elem,))

    cpdef void insert_from(self, int index, elem):
        """Insert iterable 'elem' to 'index'."""
        if not isinstance(elem, Iterable):
            raise ValueError("extend iterable only")
        cdef OrderedSet tail = self[index:]
        del self[index:]
        self.update(elem)
        self.update(tail)

    cpdef void discard(self, elem):
        """Remove element 'elem' from the 'OrderedSet' if it is present."""
        _discard(self, elem)

    cpdef void roll(self, elem, int index):
        """Roll the list to 'elem' as 'index' item."""
        if elem not in self:
            raise ValueError(f"{elem} is not in {self}")

        while elem != self[index]:
            # Rolling the list
            self.add(self.pop(0))

    cpdef object pop(self, int index = -1):
        """Remove and return the last element or an arbitrary set element.
        Raises 'KeyError' if the 'OrderedSet' is empty.
        """
        if not self:
            raise KeyError(f'{self} is empty')
        key = self[index]
        _discard(self, key)
        return key

    cpdef void remove(self, elem):
        """Remove element 'elem' from the 'set'.
        Raises KeyError if 'elem' is not contained in the set.
        """
        if elem not in self:
            raise KeyError(elem)
        _discard(self, elem)

    cpdef void clear(self):
        """Remove all elements from the 'set'."""
        cdef _Entry end = self.end
        end.next.prev = end.next = None

        # reinitialize
        self.map.clear()
        self.os_used = 0
        self.end = end = _Entry()
        end.prev = end.next = end

    cpdef OrderedSet copy(self):
        """Copy the instance."""
        return OrderedSet(self)

    def __sub__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '-' operator."""
        if not isinstance(other, Set):
            if not isinstance(other, Iterable):
                return NotImplemented
            other = OrderedSet._from_iterable(other)

        return OrderedSet._from_iterable(value for value in self if value not in other)

    def __isub__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '-=' operator."""
        if other is self:
            self.clear()
        else:
            for value in other:
                self.discard(value)
        return self

    cpdef OrderedSet intersection(self, other):
        """Implement of '&' operator."""
        return self & other

    def __and__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '&' operator."""
        if not isinstance(other, Set):
            if not isinstance(other, Iterable):
                return NotImplemented
            other = OrderedSet._from_iterable(other)

        return OrderedSet._from_iterable(value for value in self if value in other)

    def __iand__(self, it: Iterable[Any]):
        """Implement for '&=' operator."""
        for value in (self - it):
            self.discard(value)
        return self

    cpdef list ordered_intersections(self, other, bint is_loop = False):
        """Find the max length intersection with ordered detection."""
        if not isinstance(other, Iterable):
            raise TypeError("object must be iterable")

        self_list = list(self)
        other_list = list(other)
        if is_loop:
            self_list *= 2
        # Result list
        subsets = []
        cdef int matched_self
        cdef int matched_self_old = -2
        cdef int matched_other = -2
        cdef int matched_other_old = -2
        cdef OrderedSet subset = OrderedSet.__new__(OrderedSet)
        for matched_self, self_elem in enumerate(self_list):
            try:
                matched_other = other_list.index(self_elem)
            except ValueError:
                # Not found
                continue

            match_set = {matched_other_old + 1}
            if is_loop:
                match_set.add(matched_other_old + 1 - len(other_list))

            if subset and not (matched_self == matched_self_old + 1 or matched_other in match_set):
                # Add to results
                if _not_subset(subset, subsets):
                    subsets.append(subset)
                # Start new record
                subset = OrderedSet.__new__(OrderedSet)
            subset.add(self_elem)
            matched_other_old = matched_other
            matched_self_old = matched_self

        if _not_subset(subset, subsets):
            subsets.append(subset)
        return subsets

    cpdef bint isdisjoint(self, other):
        """Return true if the set has no elements in common with other.
        Sets are disjoint if and only if their intersection is the empty set.
        """
        for value in other:
            if value in self:
                return False
        return True

    cpdef bint issubset(self, other):
        """Method of '<=' operator."""
        return self <= other

    cpdef bint issuperset(self, other):
        """Method of '>=' operator."""
        return other <= self

    cpdef bint is_ordered_subset(self, other, bint is_loop = False):
        """Method of '<=' operator with ordered detection."""
        return _isorderedsubset(self, other, is_loop)

    cpdef bint is_ordered_superset(self, other, bint is_loop = False):
        """Method of '>=' operator with ordered detection."""
        return _isorderedsubset(other, self, is_loop)

    def __xor__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '^' operator."""
        if not isinstance(other, Iterable):
            return NotImplemented

        cdef OrderedSet o_other = OrderedSet._from_iterable(other)
        return (self - o_other) | (o_other - self)

    def __ixor__(self, other: Iterable[Any]):
        """Implement for '^=' operator."""
        if other is self:
            self.clear()
        else:
            if not isinstance(other, Set):
                other = self._from_iterable(other)
            for value in other:
                if value in self:
                    self.discard(value)
                else:
                    self.add(value)
        return self

    cpdef OrderedSet union(self, other):
        """Method of '|' operator."""
        return self | other

    cpdef void update(self, other):
        """Method of '|=' operator."""
        self |= other

    def __or__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '|' operator."""
        if not isinstance(other, Iterable):
            return NotImplemented
        chain = (e for s in (self, other) for e in s)
        return OrderedSet._from_iterable(chain)

    def __ior__(self, other: Iterable[Any]) -> OrderedSet:
        """Implement for '|=' operator."""
        for elem in other:
            _add(self, elem)
        return self

    ##
    # list methods
    ##
    cpdef int index(self, elem):
        """Return the index of 'elem'. Rases :class:'ValueError' if not in the OrderedSet."""
        if elem not in self:
            raise ValueError(f"{elem} is not in {type(self).__name__}")
        cdef _Entry curr = self.end.next
        cdef ssize_t index = 0
        while curr.key != elem:
            curr = curr.next
            index += 1
        return index

    cdef list _get_slice_entry(self, slice item):
        cdef ssize_t start, stop, step, slicelength
        PySlice_GetIndicesEx(item, len(self), &start, &stop, &step, &slicelength)
        result = []
        cdef ssize_t place = start
        cdef _Entry curr = self.end

        cdef ssize_t i
        if slicelength <= 0:
            pass
        elif step > 0:
            # normal forward slice
            i = 0
            while slicelength > 0:
                while i <= place:
                    curr = curr.next
                    i += 1
                result.append(curr)
                place += step
                slicelength -= 1
        else:
            # we're going backwards
            i = len(self)
            while slicelength > 0:
                while i > place:
                    curr = curr.prev
                    i -= 1
                result.append(curr)
                place += step
                slicelength -= 1
        return result

    cdef _Entry _get_index_entry(self, ssize_t index):
        cdef ssize_t _len = len(self)
        if index >= _len or (index < 0 and abs(index) > _len):
            raise IndexError("list index out of range")

        cdef _Entry curr
        if index >= 0:
            curr = self.end.next
            while index:
                curr = curr.next
                index -= 1
        else:
            index = abs(index) - 1
            curr = self.end.prev
            while index:
                curr = curr.prev
                index -= 1
        return curr

    def __getitem__(self, index):
        """Implement of 'self[index]' operator."""
        if isinstance(index, slice):
            return OrderedSet([curr.key for curr in self._get_slice_entry(index)])
        if not PyIndex_Check(index):
            raise TypeError(f"{type(self).__name__} indices must be integers, not {type(index)}")

        cdef _Entry curr = self._get_index_entry(index)
        return curr.key

    def __setitem__(self, index, value):
        """Implement of 'self[index] = value' operator."""
        cdef int i
        cdef _Entry curr
        if isinstance(index, slice):
            if not isinstance(value, Iterable):
                raise TypeError("must assign iterable to extended slice")
            value_list = list(value)
            for i, curr in enumerate(self._get_slice_entry(index)):
                curr.key = value_list[i]
            return

        if not PyIndex_Check(index):
            raise TypeError(f"{type(self).__name__} indices must be integers, not {type(index)}")
        curr = self._get_index_entry(index)
        curr.key = value

    def __delitem__(self, index):
        """Implement of 'del self[index]' operator."""
        cdef _Entry curr
        if isinstance(index, slice):
            for curr in self._get_slice_entry(index):
                self.discard(curr.key)
            return

        if not PyIndex_Check(index):
            raise TypeError(f"{type(self).__name__} indices must be integers, not {type(index)}")
        curr = self._get_index_entry(index)
        self.discard(curr.key)

    cpdef void reverse(self):
        """Reverse all elements."""
        my_iter = list(_OrderedSetReverseIterator(self))
        self.clear()
        self.update(my_iter)

    ##
    # sequence methods
    ##
    def __len__(self) -> int:
        """Implement of 'len(self)' operator."""
        return len(self.map)

    def __contains__(self, elem: Any) -> bool:
        """Implement of 'elem in self' operator."""
        return elem in self.map

    def __iter__(self) -> _OrderedSetIterator:
        """Return iterator without copy."""
        return _OrderedSetIterator.__new__(_OrderedSetIterator, self)

    def __reversed__(self) -> _OrderedSetReverseIterator:
        """Return reversed iterator without copy."""
        return _OrderedSetReverseIterator.__new__(_OrderedSetReverseIterator, self)

    def __reduce__(self):
        items = list(self)
        inst_dict = vars(self).copy()
        return type(self), (items, ), inst_dict

    def __repr__(self) -> str:
        """Implement of '!r' operator in string."""
        if not self:
            return f'{type(self).__name__}()'
        return f'{type(self).__name__}({list(self)!r})'

    def __eq__(self, other: Iterable[Any]) -> bool:
        """Implement of '==' operator."""
        if isinstance(other, Set):
            # Set is no ordered
            return set(self) == set(other)
        if isinstance(other, Iterable):
            return len(self) == len(other) and list(self) == list(other)
        return NotImplemented

    def __le__(self, other: Iterable[Any]) -> bool:
        """Implement of '<=' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) <= len(other) and set(self) <= set(other)
        elif isinstance(other, list):
            return len(self) <= len(other) and list(self) <= list(other)
        return NotImplemented

    def __lt__(self, other: Iterable[Any]) -> bool:
        """Implement of '<' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) < len(other) and set(self) < set(other)
        elif isinstance(other, list):
            return len(self) < len(other) and list(self) < list(other)
        return NotImplemented

    def __ge__(self, other: Iterable[Any]) -> bool:
        """Implement of '>=' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) >= len(other) and set(self) >= set(other)
        elif isinstance(other, list):
            return len(self) >= len(other) and list(self) >= list(other)
        return NotImplemented

    def __gt__(self, other: Iterable[Any]) -> bool:
        """Implement of '>' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) > len(other) and set(self) > set(other)
        elif isinstance(other, list):
            return len(self) > len(other) and list(self) > list(other)
        return NotImplemented
