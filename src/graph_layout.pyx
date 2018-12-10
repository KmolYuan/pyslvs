# -*- coding: utf-8 -*-
# cython: language_level=3

"""Layout functions for graph object.

The algorithms references:
+ NetworkX

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from typing import Set, Iterable
cimport cython
from cpython cimport PyDict_Contains, PyIndex_Check
from cpython.slice cimport PySlice_GetIndicesEx
from libc.math cimport M_PI, sin, cos
from graph cimport Graph


cpdef dict outer_loop_layout(Graph g, bint node_mode, double scale = 1.):
    """Layout position decided by outer loop."""
    cdef OrderedSet o_loop = _outer_loop(g)
    o_loop.roll(min(o_loop), 0)
    cdef list o_pos = _regular_polygon_pos(len(o_loop), scale)
    cdef dict pos = dict(zip(o_loop, o_pos))
    _outer_loop_layout_inner(g, o_loop, pos)

    # Last check for debug.
    if set(g.nodes) != set(pos):
        raise ValueError(
            f"the algorithm is error with {g.edges}\n"
            f"outer loop: {o_loop}\n"
            f"inner layout: {set(pos) - o_loop}\n"
            f"node {set(g.nodes) - set(pos)} are not included"
        )

    # TODO: node_mode

    return pos


cdef inline void _outer_loop_layout_inner(Graph g, OrderedSet o_loop, dict pos):
    """Layout for inner nodes of graph block."""
    # TODO: Make direction from start and end nodes.
    cdef OrderedSet nodes = set(g.nodes) - o_loop
    if not nodes:
        return

    cdef OrderedSet used_nodes = o_loop.copy()

    cdef int n, start, end
    cdef OrderedSet new_neighbors, line, inter
    while nodes:
        n = nodes.pop(0)
        if not (used_nodes & g.adj[n]):
            # Not contacted yet.
            nodes.add(n)
            continue

        line = OrderedSet.__new__(OrderedSet)
        line.add(n)
        new_neighbors = nodes & g.adj[n]
        while new_neighbors:
            # New nodes to add.
            n = new_neighbors.pop()
            line.add(n)
            new_neighbors = (nodes - line) & g.adj[n]

        # Line is ended.
        if line[0] == line[-1]:
            inter = used_nodes & g.adj[line[0]]
            start = inter.pop()
            end = inter.pop()
        else:
            inter = used_nodes & g.adj[line[0]]
            start = inter.pop()
            inter = used_nodes & g.adj[line[-1]]
            end = inter.pop()
        pos.update(zip(line, _linear_layout(pos[start], pos[end], len(line))))
        used_nodes.update(line)
        line.clear()

        nodes -= line


cdef list _regular_polygon_pos(int edge_count, double scale):
    """Return position of a regular polygon with radius 100.
    Start from bottom with _clockwise.
    """
    scale *= 5
    cdef int i
    cdef double angle = M_PI * 3 / 2
    cdef double angle_step = 2 * M_PI / edge_count
    cdef list pos = []
    for i in range(edge_count):
        pos.append((scale * cos(angle), scale * sin(angle)))
        angle -= angle_step
    return pos


cdef list _linear_layout(tuple c0, tuple c1, int count):
    """Layout position decided by equal division between two points."""
    if count < 1:
        raise ValueError(f"Invalid point number {count}")

    count += 1
    cdef double sx = (c1[0] - c0[0]) / count
    cdef double sy = (c1[1] - c0[1]) / count
    cdef int i
    return [(c0[0] + i * sx, c0[1] + i * sy) for i in range(1, count)]


cdef OrderedSet _outer_loop(Graph g):
    """Return nodes of outer loop."""
    cdef list cycles = _cycle_basis(g)
    if not cycles:
        raise ValueError(f"invalid graph has no any cycle: {g.edges}")

    cdef bint need_to_rev
    cdef int n1, n2, i, start, end
    cdef int insert_start, insert_end, replace_start, replace_end
    cdef OrderedSet c1, c2, inter
    cdef list inter_over, inter_tmp
    cdef dict inter_map
    while len(cycles) > 1:
        c1 = cycles.pop()

        for c2 in cycles:
            # Find the intersection.
            if not len(c1 & c2) >= 2:
                continue

            # Ignore subsets.
            if c1 >= c2:
                cycles.remove(c2)
                cycles.append(c1)
                break

            inter_over = c1.ordered_intersections(c2, is_loop=True)
            # Find the intersection with reversed cycle.
            c2.reverse()
            inter_tmp = c1.ordered_intersections(c2, is_loop=True)
            if len(inter_tmp) < len(inter_over):
                # Choose the longest continuous intersection.
                inter_over = inter_tmp
            inter_tmp = None

            start = -1
            end = -1
            need_to_rev = len(inter_over) % 2 == 1
            for i, inter in enumerate(inter_over):
                if not inter.is_ordered_subset(c1, is_loop=True):
                    # Intersection and cycle 1 has wrong direction.
                    inter.reverse()
                if inter.is_ordered_subset(c2, is_loop=True) == need_to_rev:
                    # Cycle 1 and cycle 2 should has different direction.
                    c2.reverse()

                # Prune cycle 2 by intersection.
                c2 -= inter[1:-1]

                # Interface nodes.
                if i == 0:
                    start = inter[0]
                end = inter[-1]

            # Roll to interface.
            c1.roll(end, -1)
            c2.roll(start, 0)

            # Insert new edges.
            insert_start = c1.index(start)
            insert_end = c1.index(end)
            replace_start = c2.index(start)
            replace_end = c2.index(end)
            if (replace_end - replace_start) > (insert_end - insert_start):
                # Cycle 2 should longer then intersection.
                del c1[insert_start:insert_end]
                c1.insert(insert_start, c2[replace_start:replace_end])

            # The cycle 2 has been merged into cycle 1.
            cycles.remove(c2)
            cycles.append(c1)
            break
        else:
            # Cycles has no contacted.
            # Find connection from edges.
            for c2 in cycles:
                inter = OrderedSet.__new__(OrderedSet)
                inter_map = {}
                for n1, n2 in g.edges:
                    if (n1 in c1 and n2 in c2) or (n1 in c2 and n2 in c1):
                        inter.add(n1)
                        inter_map[n1] = n2

                if not inter:
                    continue

                # Resort intersection.
                if not inter.is_ordered_subset(c1, is_loop=True):
                    inter = c1 & inter
                start = inter[0]
                end = inter[-1]
                insert_start = c1.index(start)
                insert_end = c1.index(end)
                replace_start = c2.index(inter_map[start])
                replace_end = c2.index(inter_map[start])
                if replace_start > replace_end:
                    c2.reverse()

                # Roll to interface.
                c1.roll(end, -1)
                c2.roll(inter_map[start], 0)

                # Merge them.
                if (replace_end - replace_start) > (insert_end - insert_start):
                    del c1[insert_start:insert_end]
                    c1.insert(insert_start, c2[replace_start:replace_end])

                # The cycle 2 has been merged into cycle 1.
                inter_map = None
                cycles.remove(c2)
                cycles.append(c1)
                break
            else:
                raise ValueError(
                    f"invalid graph: {g.edges}\n"
                    f"last one: {c1}\n"
                    f"with cycle(s): {cycles}"
                )

    return cycles.pop()


cdef inline list _cycle_basis(Graph g):
    """ Returns a list of cycles which form a basis for cycles of G.
    Reference from NetworkX.
    """
    cdef set g_nodes = set(g.nodes)
    cdef list cycles = []
    cdef int root = -1

    cdef int z, nbr, p
    cdef list stack
    cdef set zused, pn
    cdef OrderedSet cycle
    cdef dict pred, used
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
                    cycles.append([z])
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
    """Return True if 'seq1' is ordered subset of 'seq2'."""
    cdef int seq1_len = len(seq1)
    if not seq1_len <= len(seq2):
        # 'seq1' is obviously not a subset.
        return False

    cdef list seq1_list = list(seq1)
    cdef list seq2_list = list(seq2)
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
    """Return True if 'o_set' is not any subset of members."""
    cdef OrderedSet member
    for member in members:
        if member < o_set:
            # Is superset
            members.remove(member)
            return True
        if member >= o_set:
            return False
    # No subset.
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
            raise StopIteration()
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
            raise StopIteration()
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

    def __init__(self, object iterable=None):
        if iterable is None:
            return

        cdef _Entry next_e
        cdef dict map_d = self.map
        cdef _Entry end = self.end
        for elem in iterable:
            if not PyDict_Contains(map_d, elem):
                next_e = _Entry()
                next_e.key, next_e.prev, next_e.next = elem, end.prev, end
                end.prev.next = end.prev = map_d[elem] = next_e
                self.os_used += 1

    @classmethod
    def _from_iterable(cls, it) -> OrderedSet:
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
            # Rolling the list.
            self.add(self.pop(0))

    cpdef object pop(self, int index = -1):
        """Remove and return the last element or an arbitrary set element.
        Raises 'KeyError' if the 'OrderedSet' is empty.
        """
        if not self:
            raise KeyError(f'{self.__class__.__name__} is empty')
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
        return self._from_iterable(self)

    def __sub__(self, other) -> OrderedSet:
        """Implement of '-' operator."""
        if not isinstance(other, Set):
            if not isinstance(other, Iterable):
                return NotImplemented
            other = OrderedSet._from_iterable(other)

        return OrderedSet._from_iterable(value for value in self if value not in other)

    def __isub__(self, other) -> OrderedSet:
        """Implement of '-=' operator."""
        if other is self:
            self.clear()
        else:
            for value in other:
                self.discard(value)
        return self

    cpdef OrderedSet intersection(self, other):
        """Method of '&' operator."""
        return self & other

    def __and__(self, other) -> OrderedSet:
        """Implement of '&' operator."""
        if not isinstance(other, Set):
            if not isinstance(other, Iterable):
                return NotImplemented
            other = OrderedSet._from_iterable(other)

        return OrderedSet._from_iterable(value for value in self if value in other)

    def __iand__(self, it):
        """Implement of '&=' operator."""
        for value in (self - it):
            self.discard(value)
        return self

    cpdef list ordered_intersections(self, other, bint is_loop = False):
        """Find the max length intersection with ordered detection."""
        if not isinstance(other, Iterable):
            raise TypeError("object must be iterable")

        cdef list self_list = list(self)
        cdef list other_list = list(other)
        if is_loop:
            self_list *= 2

        # Result list.
        cdef list subsets = []

        cdef int matched_self
        cdef set match_set
        cdef int matched_self_old = -2
        cdef int matched_other = -2
        cdef int matched_other_old = -2
        cdef OrderedSet subset = OrderedSet.__new__(OrderedSet)
        for matched_self, self_elem in enumerate(self_list):
            try:
                matched_other = other_list.index(self_elem)
            except ValueError:
                # Not found.
                continue

            match_set = {matched_other_old + 1}
            if is_loop:
                match_set.add(matched_other_old + 1 - len(other_list))

            if subset and not (matched_self == matched_self_old + 1 or matched_other in match_set):
                # Add to results.
                if _not_subset(subset, subsets):
                    subsets.append(subset)
                # Start new record.
                subset = OrderedSet.__new__(OrderedSet)
            subset.add(self_elem)
            matched_other_old = matched_other
            matched_self_old = matched_self

        if _not_subset(subset, subsets):
            subsets.append(subset)
        return subsets

    cpdef bint isdisjoint(self, other):
        """Return True if the set has no elements in common with other.
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

    def __xor__(self, other) -> OrderedSet:
        """Implement of '^' operator."""
        if not isinstance(other, Iterable):
            return NotImplemented

        cdef OrderedSet o_other = OrderedSet._from_iterable(other)
        return (self - o_other) | (o_other - self)

    def __ixor__(self, other):
        """Implement of '^=' operator."""
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

    def __or__(self, other) -> OrderedSet:
        """Implement of '|' operator."""
        if not isinstance(other, Iterable):
            return NotImplemented
        chain = (e for s in (self, other) for e in s)
        return OrderedSet._from_iterable(chain)

    def __ior__(self, other) -> OrderedSet:
        """Implement of '|=' operator."""
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

        cdef list result = []
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
        cdef _Entry curr
        if isinstance(index, slice):
            return OrderedSet([curr.key for curr in self._get_slice_entry(index)])
        if not PyIndex_Check(index):
            raise TypeError(f"{type(self).__name__} indices must be integers, not {type(index)}")

        curr = self._get_index_entry(index)
        return curr.key

    def __setitem__(self, index, value):
        """Implement of 'self[index] = value' operator."""
        cdef int i
        cdef _Entry curr
        cdef list value_list
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
        cdef list my_iter = list(_OrderedSetReverseIterator(self))
        self.clear()
        self.update(my_iter)

    ##
    # sequence methods
    ##
    def __len__(self) -> int:
        """Implement of 'len(self)' operator."""
        return len(self.map)

    def __contains__(self, elem) -> bool:
        """Implement of 'elem in self' operator."""
        return elem in self.map

    def __iter__(self) -> _OrderedSetIterator:
        """Implement of 'iter(self)' operator. (reference method)"""
        return _OrderedSetIterator.__new__(_OrderedSetIterator, self)

    def __reversed__(self) -> _OrderedSetReverseIterator:
        """Implement of 'reversed(self)' operator."""
        return _OrderedSetReverseIterator.__new__(_OrderedSetReverseIterator, self)

    def __reduce__(self):
        items = list(self)
        inst_dict = vars(self).copy()
        return self.__class__, (items, ), inst_dict

    def __repr__(self) -> str:
        """Implement of '!r' operator in string."""
        if not self:
            return f'{self.__class__.__name__}()'
        return f'{self.__class__.__name__}({list(self)!r})'

    def __eq__(self, other) -> bool:
        """Implement of '==' operator."""
        if isinstance(other, Set):
            # Set is no ordered.
            return set(self) == set(other)
        if isinstance(other, Iterable):
            return len(self) == len(other) and list(self) == list(other)
        return NotImplemented

    def __le__(self, other) -> bool:
        """Implement of '<=' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) <= len(other) and set(self) <= set(other)
        elif isinstance(other, list):
            return len(self) <= len(other) and list(self) <= list(other)
        return NotImplemented

    def __lt__(self, other) -> bool:
        """Implement of '<' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) < len(other) and set(self) < set(other)
        elif isinstance(other, list):
            return len(self) < len(other) and list(self) < list(other)
        return NotImplemented

    def __ge__(self, other) -> bool:
        """Implement of '>=' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) >= len(other) and set(self) >= set(other)
        elif isinstance(other, list):
            return len(self) >= len(other) and list(self) >= list(other)
        return NotImplemented

    def __gt__(self, other) -> bool:
        """Implement of '>' operator."""
        if isinstance(other, (Set, OrderedSet)):
            return len(self) > len(other) and set(self) > set(other)
        elif isinstance(other, list):
            return len(self) > len(other) and list(self) > list(other)
        return NotImplemented
