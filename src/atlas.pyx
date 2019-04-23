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
from libcpp.map cimport map as cmap
from numpy cimport ndarray, int16_t
from numpy import (
    int16,
    array as np_array,
)
from graph cimport Graph
from planar_check cimport is_planar

cdef object logger = getLogger()

ctypedef unsigned int uint
ctypedef cpair[int, int] ipair
ctypedef cmap[int, int] imap


cdef int16_t[:] _labels(int16_t[:] numbers, int index, int offset, bint negative):
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


cdef inline list _picked_branch(int node, imap &limit, imap &count):
    """Return feasible node for combination."""
    cdef int pick_count = limit[node] - count[node]
    if pick_count < 1:
        return []

    # Create pool
    cdef list pool_list = []
    cdef ipair it1, it2
    for it1 in limit:
        if node == it1.first:
            continue
        if it1.second > 0:
            # Multiple links
            if count[it1.first] < it1.second:
                pool_list.append((it1.first,))
        else:
            # Contracted links
            if count[it1.first] > 0:
                continue
            for it2 in limit:
                if node == it2.first:
                    continue
                # Multiple links
                if it2.second > 0 and count[it2.first] < it2.second:
                    pool_list.append((it1.first, it2.first))

    # Check over picked
    cdef int pool_size = len(pool_list)
    if pick_count > pool_size:
        return []

    cdef int *indices = <int *>PyMem_Malloc(pick_count * sizeof(int))
    cdef int i
    for i in range(pick_count):
        indices[i] = i

    cdef bint failed = False
    cdef set types = set()
    cdef list pick_list = []
    cdef list hash_list = []
    cdef list combine_list = []

    # Combinations loop with number checking.
    cdef int n1, c1, n2, c2
    cdef tuple hash_code, hash_codes
    while True:
        # Combine
        for i in range(pick_count):
            n1 = indices[i]
            pick_list.append(pool_list[n1])
            c1 = pool_list[n1][0]
            if len(pool_list[n1]) == 1:
                # Multiple links.
                hash_code = (hash(str(c1)),)
            else:
                # Contracted links.
                c2 = pool_list[n1][1]
                hash_code = (limit[c1], hash(str(c2)))
            hash_list.append(hash_code)

        # Check if contracted link is over selected.
        if not failed and _over_count(pick_list, limit, count):
            failed = True

        # Check hash codes.
        if not failed and hash_list:
            hash_list.sort()
            hash_codes = tuple(hash_list)
            if hash_codes in types:
                failed = True
            else:
                types.add(hash_codes)

        # Collecting.
        if not failed:
            combine_list.append(tuple(pick_list))

        # Initialize.
        failed = False
        hash_list.clear()
        pick_list.clear()

        # Check combination is over.
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            PyMem_Free(indices)
            return combine_list

        # Next indexing.
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1


cdef inline int _feasible_link(imap &limit, imap &count):
    """Return next feasible multiple link, return -1 if no any matched."""
    cdef ipair it
    for it in limit:
        if it.second > 0 and count[it.first] < it.second:
            return it.first
    return -1


cdef inline bint _all_connected(imap &limit, imap &count):
    """Return True if all multiple links and contracted links is connected."""
    cdef ipair it
    for it in limit:
        if it.second < 0:
            # Contracted links.
            if count[it.first] == 0:
                return False
        else:
            # Multiple links.
            if count[it.first] != it.second:
                return False
    return True


cdef inline tuple _contracted_chain(int node, int num, set edges):
    """Get the max key and return chain."""
    cdef int m, n
    cdef int max_n = 0
    for m, n in edges:
        if m > max_n:
            max_n = m
        if n > max_n:
            max_n = n
    max_n += 1
    cdef int last_node = node
    cdef set chain = set()
    for n in range(num - 1):
        chain.add((last_node, max_n))
        last_node = max_n
        max_n += 1
    return chain, last_node


cdef inline set _dyad_patch(set edges, imap &limit):
    """Return a patched edges for contracted links."""
    cdef int last_node, u, v
    cdef set new_chain
    cdef set new_edges = edges.copy()
    cdef ipair it
    for it in limit:
        # Only for contracted links.
        if not it.second < 0 or it.second == -1:
            continue
        new_chain, last_node = _contracted_chain(it.first, abs(it.second), new_edges)
        for u, v in edges:
            # Find once.
            if it.first == u or it.first == v:
                new_edges.remove((u, v))
                if it.first == u:
                    new_edges.add((v, last_node))
                else:
                    new_edges.add((u, last_node))
                break
        new_edges.update(new_chain)
    return new_edges


# NOTE: New method
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
    set edges,
    imap &limit,
    imap *count,
    list result,
    uint no_degenerate
):
    """Test the graph."""
    # All links connected
    if not _all_connected(limit, count[0]):
        return
    # Preliminary test
    cdef Graph g = Graph.__new__(Graph, edges)
    # All connected
    if not g.is_connected():
        return
    # Cut link
    if g.has_cut_link():
        return
    # Planar graph
    if not is_planar(g):
        return

    # Result graph
    g = Graph.__new__(Graph, _dyad_patch(edges, limit))
    # Graph filter depending on degenerate option
    if no_degenerate == 0 and not g.is_degenerate():
        return
    elif no_degenerate == 1 and g.is_degenerate():
        return
    # Isomorphism
    if _is_isomorphic(g, result):
        return

    result.append(g)


cdef inline void _insert_combine(
    int node,
    tuple combine,
    set edges,
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
                edges.add((b , d))
            else:
                edges.add((d , b))
            b = d
        count[0][node] += 1
        for d in dyad:
            count[0][d] += 1


cdef void _synthesis(
    int node,
    list result,
    set edges_origin,
    imap &limit,
    imap &count_origin,
    uint no_degenerate,
    object stop_func
):
    """Recursive synthesis function."""
    # Copied edge list.
    cdef set edges
    cdef imap tmp
    cdef imap *count
    # Combinations.
    cdef int next_node
    cdef tuple combine
    cdef list branches = _picked_branch(node, limit, count_origin)
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

        _insert_combine(node, combine, edges, count)

        # Recursive or end.
        next_node = _feasible_link(limit, count[0])
        if next_node == -1:
            _test_graph(edges, limit, count, result, no_degenerate)
        else:
            _synthesis(next_node, result, edges, limit, count[0], no_degenerate, stop_func)


# NOTE: New method
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


# NOTE: New method
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


# NOTE: New method
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


# NOTE: New method
cdef inline tuple _contracted_links(tuple edges, imap &limit):
    """Combination of contracted links.
    
    pool: edges
    pick: contracted links
    
    If the edge is not picked, it represent the joint is connected directly.
    """
    cdef int pick_count = limit.size()
    if pick_count < 1:
        return ()

    # Check over picked
    cdef tuple pool_list = tuple([frozenset(edge) for edge in edges])
    cdef set single_link = set(pool_list)
    cdef int pool_size = len(pool_list)
    if pool_size - len(single_link) > pick_count > pool_size:
        return ()

    # The list including required edge(s).
    cdef object confirm_list = Counter(pool_list) - Counter(single_link)

    cdef int *indices = <int *>PyMem_Malloc(pick_count * sizeof(int))
    cdef int i
    for i in range(pick_count):
        indices[i] = i

    cdef list pick_list = []
    cdef set combine_list = set()

    # Combinations loop with number checking.
    while True:
        # Combine
        for i in range(pick_count):
            pick_list.append(pool_list[indices[i]])

        # Collecting
        if not (confirm_list - Counter(pick_list)):
            pick_list.sort()
            combine_list.add(tuple(pick_list))

        # Initialize
        pick_list.clear()

        # Check combination is over.
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            PyMem_Free(indices)
            return tuple(combine_list)

        # Next indicator
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1


# NOTE: New method
cdef inline void _graph_atlas(
    list result,
    list contracted_graph,
    imap &limit,
    uint no_degenerate,
    object stop_func
):
    """Synthesis of atlas."""
    cdef Graph g
    cdef tuple combine
    for g in contracted_graph:
        # Check if stop.
        if stop_func is not None and stop_func():
            return

        print(g)
        for combine in _contracted_links(g.edges, limit):
            print(combine)
            # TODO: contracted links combination


cdef void _splice(
    list result,
    int16_t[:] m_link,
    int16_t[:] c_link,
    uint no_degenerate,
    object stop_func
):
    """Splice multiple links by:
    
    + Connect to contracted links.
    + Connect to other multiple links.
    """
    cdef imap limit, m_limit, c_limit, count
    cdef int num
    cdef int i = 0
    for num in m_link:
        limit[i] = num
        m_limit[i] = num
        count[i] = 0
        i += 1
    for num in c_link:
        # Actual limit is 1.
        limit[i] = num
        c_limit[i] = num
        count[i] = 0
        i += 1

    # Synthesis of contracted graphs
    cdef list contracted_graphs = []
    _contracted_graph(0, contracted_graphs, [], m_limit, count, stop_func)

    # Synthesis of multiple links
    _graph_atlas(result, contracted_graphs, c_limit, no_degenerate, stop_func)

    # Origin one
    i = 0
    for num in m_link:
        count[i] = 0
        i += 1
    for num in c_link:
        count[i] = 0
        i += 1
    _synthesis(0, result, set(), limit, count, no_degenerate, stop_func)

    if len(contracted_graphs) == 0 and len(result) > 0:
        print(f"error: {len(contracted_graphs)}, {len(result)}")

    logger.debug(f"Contracted graph(s): {len(contracted_graphs)}")


cdef bint _is_isomorphic(Graph g, list result):
    """Return True if graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cdef inline list _loop_chain(int num):
    """Loop chain."""
    cdef int i
    cdef int b = 0
    cdef list chain = []
    for i in range(1, num):
        chain.append((b, i))
        b = i
    chain.append((0, b))
    return chain


cpdef tuple topo(
    object link_num_list,
    object c_j_list,
    uint no_degenerate = 1,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    link_num_list = [L2, L3, L4, ...]
    c_j_list = [NC1, NC2, NC3, ...]
    no_degenerate:
        0: only degenerate.
        1: no degenerate.
        2: all.
    stop_func: Optional[Callable[[], None]]
        stop function can check the break point and send response.
    """
    if not link_num_list:
        return [], 0.

    # NumPy array type.
    cdef ndarray[int16_t, ndim=1] link_num = np_array(link_num_list, ndmin=1, dtype=int16)
    logger.debug(f"Assortments: {link_num} {c_j_list}")

    # Initial time.
    cdef double t0 = time()

    # Multiple links.
    cdef int16_t[:] m_link = _labels(link_num, 3, 1, False)

    # Synthesis of contracted link and multiple link combination.
    cdef int16_t[:] c_j = np_array(c_j_list, ndmin=1, dtype=int16)

    cdef Graph g
    cdef list result = []
    if len(m_link) == 0:
        # Single loop (Special case).
        result.append(Graph.__new__(Graph, _loop_chain(link_num[0])))
    else:
        _splice(result, m_link, _labels(c_j, 1, 0, True), no_degenerate, stop_func)

    # Return graph list and time.
    logger.debug(f"Count: {len(result)}")
    return result, (time() - t0)
