# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from itertools import product, combinations
from time import time
from cpython cimport bool
from libcpp.map cimport map
from numpy cimport ndarray, int64_t
from numpy import array as np_array
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

    return np_array(cj_list, dtype=int)


cdef ndarray[int64_t, ndim=1] labels(ndarray[int64_t, ndim=1] numbers, int index, int offset):
    """Generate labels from numbers."""
    cdef int i, num
    cdef list labels = []
    for num in numbers[offset:]:
        for i in range(num):
            labels.append(index)
        index += 1
    return np_array(labels, dtype=int)


cdef bool is_isomorphic(Graph g, list result):
    """Return True is graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cdef void splice(
    list result,
    ndarray[int64_t, ndim=1] m_link,
    ndarray[int64_t, ndim=1] c_link,
    object stop_func = None
) except *:
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

    # Number of all connection.
    cdef int roads = <int>((sum(m_link) + 2 * len(c_link)) / 2)
    print(roads)

    # TODO: Synthesis of multiple links.
    # TODO: Break point of stop_func.
    cdef int ml, p, b, num3, pick_count
    cdef tuple combine, dyad, edge
    cdef list pool = []
    cdef set edges = set()
    for ml, num1 in limit:
        pool.clear()
        for p, num2 in limit:
            if p == ml:
                continue
            if num2 < 0:
                if count[p] > 0:
                    continue
                # If p is a contracted link.
                for b, num3 in limit:
                    if num3 < 0 or num3 - count[b] <= 0:
                        continue
                    # If b is a multiple link.
                    pool.append((p, b))
            else:
                if num2 - count[p] <= 0:
                    continue
                # If p is a multiple link.
                pool.append((p,))

        pick_count = num1 - count[ml]
        if pick_count <= 0:
            continue
        for combine in combinations(pool, num1 - count[ml]):
            # TODO: Add connection.
            for dyad in combine:
                b = ml
                for p in dyad:
                    if b < p:
                        edge = (b, p)
                    else:
                        edge = (p, b)
                    edges.add(edge)
                    count[b] += 1
                    count[p] += 1
                    b = p


cpdef tuple topo(
    object link_num_,
    bool degenerate = True,
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
    cdef ndarray[int64_t, ndim=1] link_num = np_array(link_num_, dtype=int)

    # Multiple links.
    cdef ndarray[int64_t, ndim=1] m_link = labels(link_num, 3, 1)

    # Synthesis of contracted link and multiple link combination.
    cdef ndarray[int64_t, ndim=1] c_j
    cdef list result = []
    for c_j in contracted_link(link_num):
        # TODO: Limitation of job_func.
        # job_func(str(c_j), len(m_link) * len(c_j))
        splice(result, m_link, -labels(c_j, 1, 0), stop_func)

    cdef Graph g
    print([g.edges for g in result])
    print(f"Count: {len(result)}")
    print(f"Time: {time() - t0:.04f}")
    # Return graph list and time.
    return result, (time() - t0)
