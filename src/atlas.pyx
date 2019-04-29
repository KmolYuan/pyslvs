# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True, cdivision=True

"""Structure _synthesis.

The algorithms references:
+ On the Number Synthesis of Kinematic Chains
    + author: Hong-Sen Yan, Yu-Ting Chiu
    + Mechanism and Machine Theory, Volume 89, 2015, Pages 128-144, ISSN 0094-114X
    + https://doi.org/10.1016/j.mechmachtheory.2014.08.012

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

from time import time
from logging import getLogger
from collections import Counter
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libcpp.pair cimport pair as cpair
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from numpy cimport int16_t
from numpy import (
    int16,
    array as np_array,
    abs as np_abs,
)
from graph cimport Graph, link_assortments
from planar_check cimport is_planar

ctypedef unsigned int uint
ctypedef cpair[int, int] ipair
ctypedef cmap[int, int] imap

cdef object logger = getLogger()


cdef inline int16_t[:] _labels(
    int16_t[:] numbers,
    int index,
    int offset,
    bint negative
):
    """Generate labels from numbers."""
    cdef int i, num
    cdef list labels = []
    for num in numbers[offset:]:
        for i in range(num):
            if negative:
                labels.append(-index)
            else:
                labels.append(index)
        index += 1
    return np_array(labels, dtype=int16)


cdef inline bint _over_count(list pick_list, imap &limit, imap &count):
    """Return True if it is a feasible pick list."""
    cdef int n
    cdef tuple candidate
    cdef imap pre_count
    for candidate in pick_list:
        for n in candidate:
            pre_count[n] += 1
            if limit[n] < 0:
                # Contracted links.
                if pre_count[n] + count[n] > 1:
                    return True
            else:
                # Multiple links.
                if pre_count[n] + count[n] > limit[n]:
                    return True
    return False


cdef inline int _feasible_link(imap &limit, imap &count):
    """Return next feasible multiple link, return -1 if no any matched."""
    cdef ipair it
    for it in limit:
        if it.second > 0 and count[it.first] < it.second:
            return it.first
    return -1


cdef inline bint _is_isomorphic(Graph g, list result):
    """Return True if graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cdef inline void _test_contracted_graph(
    Graph g,
    imap &limit,
    imap *count,
    list result
):
    """Test the contracted graph."""
    # All connected
    if not g.is_connected():
        return
    # Cut link
    if g.has_cut_link():
        return
    # Planar graph
    if not is_planar(g):
        return
    # Isomorphism
    if _is_isomorphic(g, result):
        return

    result.append(g)


cdef inline void _test_graph(
    Graph g,
    list result,
    uint no_degenerate
):
    """Test result graph."""
    # Graph filter depending on degenerate option
    if no_degenerate == 0 and not g.is_degenerate():
        return
    elif no_degenerate == 1 and g.is_degenerate():
        return
    # Isomorphism
    if _is_isomorphic(g, result):
        return
    result.append(g)


cdef inline list _picked_multi_branch(int node, imap &limit, imap &count):
    """Return feasible node for contracted graph combination."""
    cdef int pick_count = limit[node] - count[node]
    if pick_count < 1:
        return []

    # Create pool
    cdef list pool_list = []
    cdef int i, margin
    cdef ipair it1, it2
    for it1 in limit:
        if node == it1.first:
            continue
        if it1.second <= 0:
            continue
        # Multiple links
        margin = it1.second - count[it1.first]
        if margin <= 0:
            continue
        for i in range(margin):
            pool_list.append((it1.first,))

    # Check over picked
    cdef int pool_size = len(pool_list)
    if pick_count > pool_size:
        return []

    cdef int *indices = <int *>PyMem_Malloc(pick_count * sizeof(int))
    for i in range(pick_count):
        indices[i] = i

    cdef list pick_list = []
    cdef list combine_list = []

    # Combinations loop with number checking.
    cdef int n1, n2
    while True:
        # Combine
        for i in range(pick_count):
            pick_list.append(pool_list[indices[i]])

        # Check if contracted link is over selected.
        if not _over_count(pick_list, limit, count):
            combine_list.append(tuple(pick_list))

        # Initialize
        pick_list.clear()

        # Check combination is over.
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            PyMem_Free(indices)
            return combine_list

        # Next indicator
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1


cdef inline void _insert_edges(
    int node,
    tuple combine,
    list edges,
    imap *count
):
    """Insert combinations."""
    # Collecting to edges.
    cdef int b, d
    cdef tuple dyad
    for dyad in combine:
        b = node
        for d in dyad:
            if b < d:
                edges.append((b , d))
            else:
                edges.append((d , b))
            b = d
        count[0][node] += 1
        for d in dyad:
            count[0][d] += 1


cdef void _contracted_graph(
    int node,
    list result,
    list edges_origin,
    imap &limit,
    imap &count_origin,
    object stop_func
):
    """Synthesis of contracted graphs."""
    # Copied edge list.
    cdef list edges
    cdef imap tmp
    cdef imap *count
    # Combinations.
    cdef int next_node
    cdef tuple combine
    cdef list branches = _picked_multi_branch(node, limit, count_origin)
    cdef bint multi_case = len(branches) > 1
    for combine in branches:
        # Check if stop.
        if stop_func is not None and stop_func():
            return

        if multi_case:
            edges = edges_origin.copy()
            tmp = count_origin
            count = &tmp
        else:
            edges = edges_origin
            count = &count_origin

        _insert_edges(node, combine, edges, count)

        # Recursive or end.
        next_node = _feasible_link(limit, count[0])
        if next_node == -1:
            _test_contracted_graph(Graph.__new__(Graph, edges), limit, count, result)
        else:
            _contracted_graph(next_node, result, edges, limit, count[0], stop_func)


cdef inline void _dyad_insert(Graph g, frozenset edge, int amount):
    """Insert dyad to the graph."""
    if amount < 1:
        return

    cdef int last_num = max(g.nodes) + 1
    cdef int n1, n2
    n1, n2 = edge

    cdef list path = [n1]
    path.extend(range(last_num, last_num + amount))
    path.append(n2)

    g.remove_edge(n1, n2)
    g.add_path(path)


cdef inline void _permute_combine(
    int16_t[:] limit,
    list combine_list,
    tuple pick_list
):
    """Permutation of combined list."""
    cdef int n = len(limit)
    if n < 1:
        return

    cdef vector[int] indices = range(n)
    cdef int *cycles = <int *>PyMem_Malloc(n * sizeof(int))
    cdef int16_t[:] pool = np_abs(limit)

    cdef set permute_list = set()

    cdef int i, j
    for i, j in enumerate(range(n, 0, -1)):
        cycles[i] = j

    permute_list.add(tuple(pool[indices[i]] for i in range(n)))

    cdef vector[int].iterator it2 = indices.begin()

    cdef int tmp
    while True:
        for i in reversed(range(n)):
            cycles[i] -= 1
            if cycles[i] == 0:
                tmp = indices[i]
                indices.erase(it2 + i)
                indices.push_back(tmp)
                cycles[i] = n - i
            else:
                j = cycles[i]
                tmp = indices[i]
                indices[i] = indices[n - j]
                indices[n - j] = tmp
                permute_list.add(tuple(pool[indices[i]] for i in range(n)))
                break
        else:
            break

    PyMem_Free(cycles)

    cdef tuple tmp_array
    for tmp_array in permute_list:
        combine_list.append(tuple(zip(pick_list, tmp_array)))


cdef inline list _contracted_links(tuple edges, int16_t[:] limit):
    """Combination of contracted links.
    
    pool: edges
    pick: contracted links
    
    If the edge is not picked, it represent the joint is connected directly.
    """
    # Number of contracted links
    cdef int pick_count = len(limit)
    if pick_count < 1:
        return []

    # Check over picked
    cdef tuple pool_list = tuple([frozenset(edge) for edge in edges])
    cdef object pool = Counter(pool_list)
    # The list including required edge(s).
    cdef object confirm_list = pool - Counter(set(pool_list))
    # Simplified the pool
    pool_list = tuple((pool - confirm_list).elements())
    cdef int confirm_size = sum(confirm_list.values())
    if confirm_size > pick_count or pick_count > len(edges):
        return []

    pick_count -= confirm_size
    cdef int pool_size = len(pool_list)
    cdef int *indices = <int *>PyMem_Malloc(pick_count * sizeof(int))
    cdef int i
    for i in range(pick_count):
        indices[i] = i

    cdef object pick_list = Counter()
    cdef set combine_set = set()

    # Combinations loop with number checking.
    while True:
        # Combine
        for i in range(pick_count):
            pick_list[pool_list[indices[i]]] += 1

        # Collecting
        combine_set.add(tuple(pick_list.elements()))

        # Initialize
        pick_list.clear()

        # Check combination is over.
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            break

        # Next indicator
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1

    PyMem_Free(indices)

    cdef list combine_list = []
    pool_list = tuple(confirm_list.elements())
    cdef tuple picked
    for picked in combine_set:
        _permute_combine(limit, combine_list, pool_list + picked)
    return combine_list


cdef inline void _graph_atlas(
    list result,
    list contracted_graph,
    int16_t[:] limit,
    uint no_degenerate,
    object stop_func
):
    """Synthesis of atlas."""
    cdef int n
    cdef Graph cg, g
    cdef tuple combine
    cdef frozenset edge
    for cg in contracted_graph:
        # Check if stop.
        if stop_func is not None and stop_func():
            return

        for combine in _contracted_links(cg.edges, limit):
            g = Graph.__new__(Graph, cg.edges)
            for edge, n in combine:
                _dyad_insert(g, edge, n)
            _test_graph(g, result, no_degenerate)


cdef inline list _loop_chain(int num):
    """Loop chain of num."""
    cdef int i
    cdef int b = 0
    cdef list chain = []
    for i in range(1, num):
        chain.append((b, i))
        b = i
    chain.append((0, b))
    return chain


cpdef list contracted_graph(object link_num_list, object stop_func = None):
    """Get contracted graph by link assortment."""
    if not link_num_list:
        return []

    # Initial time
    cdef double t0 = time()
    cdef int16_t[:] link_num = np_array(link_num_list, ndmin=1, dtype=int16)
    logger.debug(f"Link assortment: {list(link_num)}")

    # Multiple links
    cdef int16_t[:] m_link = _labels(link_num, 3, 1, False)

    cdef imap m_limit, count
    cdef int num
    cdef int i = 0
    for num in m_link:
        m_limit[i] = num
        count[i] = 0
        i += 1

    # Synthesis of contracted graphs
    cdef list cg_list = []
    _contracted_graph(0, cg_list, [], m_limit, count, stop_func)

    logger.debug(f"Contracted graph(s): {len(cg_list)}, time: {time() - t0}")
    return cg_list


cpdef list topo(
    list cg_list,
    object c_j_list,
    uint no_degenerate = 1,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    cg_list: Contracted graph list (List[Graph]).
    c_j_list = [NC1, NC2, NC3, ...]
    no_degenerate:
        0: only degenerate.
        1: no degenerate.
        2: all.
    stop_func: Optional[Callable[[], None]]
        stop function can check the break point and send response.
    """
    # Initial time
    cdef double t0 = time()
    logger.debug(f"Contracted link assortment: {list(c_j_list)}")

    # Synthesis of contracted link and multiple link combination.
    cdef int16_t[:] c_j = np_array(c_j_list, ndmin=1, dtype=int16)

    cdef list result = []
    if not cg_list:
        if 1 not in c_j_list:
            return []
        # Single loop - ring graph (special case)
        result.append(Graph.__new__(Graph, _loop_chain(c_j_list.index(1))))
        logger.debug(f"Count: {len(result)}")
        return result

    # Multiple links
    cdef int16_t[:] m_link = np_array(link_assortments(cg_list[0]), ndmin=1, dtype=int16)
    m_link = _labels(m_link, 3, 1, False)

    cdef imap m_limit, count
    cdef int num
    cdef int i = 0
    for num in m_link:
        m_limit[i] = num
        count[i] = 0
        i += 1

    # Synthesis of multiple links
    _graph_atlas(result, cg_list, _labels(c_j, 1, 0, True), no_degenerate, stop_func)

    # Return graph list and time
    logger.debug(f"Count: {len(result)}, time: {time() - t0}")
    return result
