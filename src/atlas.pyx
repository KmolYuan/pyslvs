# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis.

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
    # Contracted link.
    cdef int n_c_max, n_c_min
    n_c_max, n_c_min = n_c(link_num)
    if n_c_max < n_c_min:
        n_c_max, n_c_min = n_c_min, n_c_max

    # NL2 - Nc + 2
    cdef int i_max = link_num[0] - n_c_min + 2
    cdef int i_min = link_num[0] - n_c_max + 2

    cdef tuple m
    cdef int count, factor, index
    cdef list cj_list = []
    for m in product(range(i_max + 1), repeat=i_max - 1):
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

    return np_array(cj_list, dtype=int64)


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
    map[int, int] limit,
    map[int, int] count,
    list pool_list,
    int pick_count
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
        for n2 in pair:
            if n2 == n1:
                continue
            if count[n2] == 0 and limit[n2] == limit[n1]:
                if counter > pick_count:
                    return True
                counter += 1
    return False


cdef list pool(int node, map[int, int] limit, map[int, int] count):
    """Return feasible node for combination."""
    cdef int pick_count = limit[node] - count[node]
    if pick_count <= 0:
        return []

    # Create pool.
    cdef int n1, n2, c1, c2
    cdef list pool_list = []
    for n1, c1 in limit:
        if node == n1 or is_same_pool(n1, limit, count, pool_list, pick_count):
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
                if c2 > 0 and count[n2] == 0:
                    pool_list.append((n1, n2))

    # Combine.
    # TODO: Need to avoid the point that has picked.
    cdef int pool_size = len(pool_list)
    if pick_count > pool_size:
        return []
    cdef list combine_list = []
    cdef vector[int] indices = range(pick_count)
    combine_list.append(tuple(pool_list[n1] for n1 in indices))
    while True:
        for n1 in reversed(range(pick_count)):
            if indices[n1] != n1 + pool_size - pick_count:
                break
        else:
            return combine_list
        indices[n1] += 1
        for n2 in range(n1 + 1, pick_count):
            indices[n2] = indices[n2 - 1] + 1
        combine_list.append(tuple(pool_list[n1] for n1 in indices))


cdef inline int feasible_link(map[int, int] limit, map[int, int] count):
    """Return next feasible multiple link, return -1 if no any matched."""
    cdef int n, c
    for n, c in limit:
        if c > 0 and count[n] < c:
            return n
    return -1


cdef inline bool all_connected(set edges, map[int, int] limit, map[int, int] count):
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


cdef inline set dyad_patch(set edges, map[int, int] limit):
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
    map[int, int] limit,
    map[int, int] count_origin
):
    """Recursive synthesis function."""
    # Copied edge list.
    cdef set edges
    # Combinations.
    cdef int b, d, next_node
    cdef Graph g
    cdef tuple combine, dyad
    for combine in pool(node, limit, count_origin):
        edges = edges_origin.copy()
        count = count_origin
        # Collecting to edges.
        for dyad in combine:
            b = node
            for d in dyad:
                if b < d:
                    edges.add((b , d))
                else:
                    edges.add((d , b))
                count[b] += 1
                count[d] += 1
                b = d
        # Recursive or end.
        next_node = feasible_link(limit, count)
        if next_node == -1:
            # Is all links connected.
            if not all_connected(edges, limit, count):
                continue
            # Transform function of contracted links.
            g = Graph(dyad_patch(edges, limit))
            # Is graph all connected.
            if not g.is_connected():
                continue
            # Collecting to result.
            result.append(g)
        else:
            synthesis(next_node, result, edges, limit, count)


cdef void splice(
    list result,
    ndarray[int64_t, ndim=1] m_link,
    ndarray[int64_t, ndim=1] c_link
):
    """Splice multiple links by:
    
    + Connect to contracted links.
    + Connect to other multiple links.
    """
    cdef map[int, int] limit, count
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
    synthesis(0, result, edges, limit, count)


cdef bool is_isomorphic(Graph g, list result):
    """Return True is graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cpdef tuple topo(
    object link_num_,
    bool no_degenerate = True,
    object job_func = None,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    link_num_ = [L2, L3, L4, ...]
    links = [[number_code]: joint_number, ...]
    """
    if not link_num_:
        return [], 0.

    # Initial time.
    cdef double t0 = time()

    # NumPy array type.
    cdef ndarray[int64_t, ndim=1] link_num = np_array(link_num_, dtype=int64)

    # Multiple links.
    cdef ndarray[int64_t, ndim=1] m_link = labels(link_num, 3, 1)

    # Start job.
    cdef ndarray[int64_t, ndim=2] c_links = contracted_link(link_num)
    job_func(str(c_links), len(c_links))

    # Synthesis of contracted link and multiple link combination.
    cdef ndarray[int64_t, ndim=1] c_j
    cdef list result = []
    for c_j in c_links:
        if stop_func and stop_func():
            break
        splice(result, m_link, -labels(c_j, 1, 0))

    cdef Graph g
    cdef list result_no_repeat = []
    for g in result:
        # If has triangles.
        if g.has_triangles() and no_degenerate:
            continue
        if is_isomorphic(g, result_no_repeat):
            continue
        result_no_repeat.append(g)

    print(f"Count: {len(result_no_repeat)}")
    print(f"Time: {time() - t0:.04f}")
    # Return graph list and time.
    return result_no_repeat, (time() - t0)
