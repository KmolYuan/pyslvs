# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis.

The algorithms references:
+ On the Number Synthesis of Kinematic Chains
    + author: Hong-Sen Yan, Yu-Ting Chiu
    + Mechanism and Machine Theory, Volume 89, 2015, Pages 128-144, ISSN 0094-114X
    + https://doi.org/10.1016/j.mechmachtheory.2014.08.012

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from itertools import product
from time import time
from cpython cimport bool
from libcpp.vector cimport vector
from libcpp.map cimport map
from numpy cimport ndarray, int64_t
from numpy import (
    array as np_array,
    int64,
)
from graph cimport Graph

ctypedef map[int, int] map_int


cdef int j_m(ndarray[int64_t, ndim=1] link_num):
    """Return value of Jm."""
    cdef int num
    cdef int i = 3
    cdef double c = 0
    for num in link_num[1:]:
        c += i / 2 * num
        i += 1
    return <int>c


cdef int jp_m(ndarray[int64_t, ndim=1] link_num):
    """Return value of J'm."""
    # Number of multiple links.
    cdef int n_m = sum(link_num[1:])
    if n_m <= 1:
        return 0
    elif n_m % 2 == 0:
        return <int>((3 * (n_m - 1) - 1) / 2)
    else:
        return <int>((3 * (n_m - 1) - 2) / 2)


cdef tuple n_c(ndarray[int64_t, ndim=1] link_num):
    """Return all values of Nc."""
    cdef int j_m_v = j_m(link_num)
    cdef int jp_m_v = jp_m(link_num)
    return max(1, j_m_v - jp_m_v), min(link_num[0], j_m_v)


cdef ndarray[int64_t, ndim=2] contracted_link(ndarray[int64_t, ndim=1] link_num):
    """Generate the contracted link assortments."""
    print("Link assortment:", link_num)
    # Contracted link.
    cdef int n_c_max, n_c_min
    n_c_max, n_c_min = n_c(link_num)
    if n_c_max < n_c_min:
        n_c_max, n_c_min = n_c_min, n_c_max

    # NL2 - Nc + 2
    cdef int i_max = link_num[0] - n_c_min + 2
    print("i max:", f"{link_num[0]} - {n_c_min} + 2 =", i_max)

    cdef int count, factor, index
    cdef tuple m
    cdef list cj_list = []
    for m in product(range(link_num[0] + 1), repeat=i_max - 1):
        # First formula.
        if not (n_c_min <= sum(m) <= n_c_max):
            continue
        # Second formula.
        count = 0
        index = 1
        for factor in m:
            count += index * factor
            index += 1
        if count == link_num[0]:
            cj_list.append(m)

    return np_array(cj_list, ndmin=2, dtype=int64)


cdef ndarray[int64_t, ndim=1] labels(
    ndarray[int64_t, ndim=1] numbers,
    int index,
    int offset
):
    """Generate labels from numbers."""
    cdef int i, num
    cdef list labels = []
    for num in numbers[offset:]:
        for i in range(num):
            labels.append(index)
        index += 1
    return np_array(labels, dtype=int64)


cdef inline bool is_same_pool(
    int n1,
    map_int &limit,
    map_int &count,
    list pool_list,
    int pick_count,
    int index
):
    """Return True if the multiple link is a duplicate status."""
    # self status should be unconnected.
    if count[n1] != 0:
        return False
    # Another status.
    cdef int n2
    cdef tuple pair
    cdef int counter = 0
    for pair in pool_list:
        if not index < len(pair):
            continue
        n2 = pair[index]
        # Same type unlink nodes.
        if limit[n2] != limit[n1] or n2 == n1 or count[n2] != 0:
            continue
        if counter > pick_count:
            return True
        counter += 1
    return False


cdef inline bool feasible_pick_list(list pick_list, map_int &limit, map_int &count):
    """Return True if it is a feasible pick list."""
    cdef int n
    cdef tuple combine
    cdef map_int pre_count
    for combine in pick_list:
        for n in combine:
            if pre_count.find(n) != pre_count.end():
                pre_count[n] = 1
            else:
                pre_count[n] += 1
            if limit[n] < 0:
                # Contracted links.
                if pre_count[n] + count[n] > 1:
                    return False
            else:
                # Multiple links.
                if pre_count[n] + count[n] > limit[n]:
                    return False
    return True


cdef inline list picked_branch(int node, map_int &limit, map_int &count):
    """Return feasible node for combination."""
    # TODO: Need to be optimized.
    cdef int pick_count = limit[node] - count[node]
    if pick_count <= 0:
        return []

    # Create pool.
    cdef int n1, n2, c1, c2
    cdef list pool_list = []
    for n1, c1 in limit:
        if node == n1:
            continue
        if is_same_pool(n1, limit, count, pool_list, pick_count, 0):
            continue
        if c1 > 0:
            # Multiple links.
            if count[n1] < c1:
                pool_list.append((n1,))
        else:
            # Contracted links.
            if count[n1] > 0:
                continue
            for n2, c2 in limit:
                if node == n2:
                    continue
                if is_same_pool(n2, limit, count, pool_list, pick_count, 1):
                    continue
                # Multiple links.
                if c2 > 0 and count[n2] < c2:
                    pool_list.append((n1, n2))

    # Over picked (error graph).
    cdef int pool_size = len(pool_list)
    if pick_count > pool_size:
        return []

    # Combine.
    cdef tuple p
    cdef list combine_list = []
    cdef vector[int] indices = range(pick_count)

    # Combinations loop with number checking.
    cdef list pick_list = [pool_list[n1] for n1 in indices]
    # Check if contracted link is over selected.
    if feasible_pick_list(pick_list, limit, count):
        combine_list.append(tuple(pick_list))
    while True:
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            return combine_list
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1
        pick_list.clear()
        for n1 in indices:
            pick_list.append(pool_list[n1])
        # Check if contracted link is over selected.
        if not feasible_pick_list(pick_list, limit, count):
            continue
        combine_list.append(tuple(pick_list))


cdef inline int feasible_link(map_int &limit, map_int &count):
    """Return next feasible multiple link, return -1 if no any matched."""
    cdef int n, c
    for n, c in limit:
        if c > 0 and count[n] < c:
            return n
    return -1


cdef inline bool all_connected(map_int &limit, map_int &count):
    """Return True if all multiple links and contracted links is connected."""
    cdef int n, c
    for n, c in limit:
        if c < 0:
            # Contracted links.
            if count[n] == 0:
                return False
        else:
            # Multiple links.
            if count[n] != c:
                return False
    return True


cdef inline tuple contracted_chain(int node, int num, set edges):
    """Get the max key and return chain."""
    cdef int m, n
    cdef int max_n = 0
    for m, n in edges:
        if m > max_n:
            max_n = m
        if n > max_n:
            max_n = n
    max_n += 1
    cdef int b = node
    cdef set chain = set()
    for n in range(num - 1):
        chain.add((b, max_n))
        b = max_n
        max_n += 1
    return chain, b


cdef inline set dyad_patch(set edges, map_int &limit):
    """Return a patched edges for contracted links."""
    cdef int n, c, b, u, v
    cdef set new_chain
    cdef set new_edges = edges.copy()
    for n, c in limit:
        # Only for contracted links.
        if not c < 0 or c == -1:
            continue
        new_chain, b = contracted_chain(n, abs(c), new_edges)
        for u, v in edges:
            # Find once.
            if n == u or n == v:
                new_edges.remove((u, v))
                if n == u:
                    new_edges.add((v, b))
                else:
                    new_edges.add((u, b))
                break
        new_edges.update(new_chain)
    return new_edges


cdef void synthesis(
    int node,
    list result,
    set edges_origin,
    map_int &limit,
    map_int &count_origin,
    object stop_func
):
    """Recursive synthesis function."""
    # Copied edge list.
    cdef set edges
    cdef map_int tmp
    cdef map_int *count
    # Combinations.
    cdef int b, d, next_node
    cdef Graph g
    cdef tuple combine
    cdef tuple dyad
    cdef list branches = picked_branch(node, limit, count_origin)
    for combine in branches:
        # Check if stop.
        if stop_func and stop_func():
            return
        if len(branches) > 1:
            edges = edges_origin.copy()
            tmp = count_origin
            count = &tmp
        else:
            edges = edges_origin
            count = &count_origin
        # Collecting to edges.
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
        # Recursive or end.
        next_node = feasible_link(limit, count[0])
        if next_node == -1:
            # Is all links connected.
            if not all_connected(limit, count[0]):
                continue
            # Transform function of contracted links.
            g = Graph(dyad_patch(edges, limit))
            # Is graph all connected.
            if not g.is_connected():
                continue
            # Is graph has cut link.
            if g.has_cut_link():
                continue
            # Collecting to result.
            result.append(g)
        else:
            synthesis(next_node, result, edges, limit, count[0], stop_func)


cdef void splice(
    list result,
    ndarray[int64_t, ndim=1] m_link,
    ndarray[int64_t, ndim=1] c_link,
    object stop_func,
):
    """Splice multiple links by:
    
    + Connect to contracted links.
    + Connect to other multiple links.
    """
    cdef map_int limit, count
    cdef int num1, num2
    cdef int i = 0
    for num1 in m_link:
        limit[i] = num1
        count[i] = 0
        i += 1
    for num2 in c_link:
        # Actual limit is 1.
        limit[i] = num2
        count[i] = 0
        i += 1

    # Synthesis of multiple links.
    cdef set edges = set()
    synthesis(0, result, edges, limit, count, stop_func)


cdef bool is_isomorphic(Graph g, list result):
    """Return True if graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cdef inline list loop_chain(int num):
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
    object link_num_,
    bool no_degenerate = True,
    object job_func = None,
    object step_func = None,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    link_num_ = [L2, L3, L4, ...]
    job_func: Optional[Callable[[List[int]], None]]
    step_func: Optional[Callable[[], None]]
    stop_func: Optional[Callable[[], None]]
    """
    if not link_num_:
        return [], 0.

    # NumPy array type.
    cdef ndarray[int64_t, ndim=1] link_num = np_array(link_num_, dtype=int64)

    # Single loop (Special case).
    if len(link_num) == 1:
        return [Graph(loop_chain(link_num[0]))], 0.

    # Initial time.
    cdef double t0 = time()

    # Multiple links.
    cdef ndarray[int64_t, ndim=1] m_link = labels(link_num, 3, 1)

    # Start job.
    cdef ndarray[int64_t, ndim=2] c_links = contracted_link(link_num)

    # Synthesis of contracted link and multiple link combination.
    cdef ndarray[int64_t, ndim=1] c_j
    if job_func:
        for c_j in c_links:
            job_func(list(c_j))

    cdef list result = []
    for c_j in c_links:
        print(c_j)
        splice(result, m_link, -labels(c_j, 1, 0), stop_func)
        if step_func:
            step_func()

    print(f"Done. Start compare results ({len(result)}) ...")

    cdef Graph g
    cdef list result_no_repeat = []
    for g in result:
        if stop_func and stop_func():
            break
        # If graph is degenerate.
        if g.is_degenerate() and no_degenerate:
            continue
        # If graph is repeated.
        if is_isomorphic(g, result_no_repeat):
            continue
        result_no_repeat.append(g)
    if step_func:
        step_func()

    print(f"Count: {len(result_no_repeat)}")
    print(f"Time: {time() - t0:.04f}")
    # Return graph list and time.
    return result_no_repeat, (time() - t0)
