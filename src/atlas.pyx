# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from itertools import product, permutations
from time import time
from cpython cimport bool
from numpy cimport ndarray
from numpy import array as np_array
from graph cimport Graph


cdef int j_m(ndarray[int, ndim=1] link_num):
    """Return value of Jm."""
    cdef int num
    cdef int i = 3
    cdef double c = 0
    for num in link_num[1:]:
        c += i / 2 * num
        i += 1
    return <int>c


cdef int jp_m(ndarray[int, ndim=1] link_num):
    """Return value of J'm."""
    # Number of multiple links.
    cdef int n_m = sum(link_num[1:])
    if n_m <= 1:
        return 0
    elif n_m % 2 == 0:
        return <int>((3 * (n_m - 1) - 1) / 2)
    else:
        return <int>((3 * (n_m - 1) - 2) / 2)


cdef tuple n_c(ndarray[int, ndim=1] link_num):
    """Return all values of Nc."""
    cdef int j_m_v = j_m(link_num)
    cdef int jp_m_v = jp_m(link_num)
    return max(1, j_m_v - jp_m_v), min(link_num[0], j_m_v)


cdef ndarray[int, ndim=2] contracted_link(ndarray[int, ndim=1] link_num):
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

    return np_array(cj_list)


cdef ndarray[int, ndim=1] labels(ndarray[int, ndim=1] numbers, int index, int offset):
    """Generate labels from numbers."""
    cdef int i, num
    cdef list labels = []
    for num in numbers[offset:]:
        for i in range(num):
            labels.append(index)
        index += 1
    return np_array(labels)


cdef void splice(list result, ndarray[int, ndim=1] m_link, ndarray[int, ndim=1] c_link):
    """Splice multiple links by:
    
    + Connect to contracted links.
    + Connect to other multiple links.
    """
    print("Contracted links:", c_link)
    cdef dict labels = {}
    cdef int num
    cdef int i = 0
    for num in m_link:
        labels[i] = num
        i += 1
    for num in c_link:
        labels[i] = num
        i += 1
    print(labels)


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
    cdef ndarray[int, ndim=1] link_num = np_array(link_num_, dtype=int)

    # Multiple links.
    cdef ndarray[int, ndim=1] m_link = labels(link_num, 3, 1)
    print("Multiple link:", m_link)

    # Synthesis of contracted link and multiple link combination.
    cdef int ml
    cdef ndarray[int, ndim=1] c_j
    cdef list result = []
    for c_j in contracted_link(link_num):
        splice(result, m_link, -labels(c_j, 1, 0))

    print(f"Time: {time() - t0:.04f}")
    # Return graph list and time.
    return result, (time() - t0)
