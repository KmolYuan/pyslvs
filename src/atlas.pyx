# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from itertools import combinations, product
from time import time
from cpython cimport bool
from numpy cimport ndarray
from numpy import array as np_array
from graph cimport Graph


cdef inline bool is_isomorphic(Graph graph1, list answer):
    """Return True if the graph is isomorphic with list."""
    cdef Graph graph2
    for i, graph2 in enumerate(answer):
        if graph1.is_isomorphic(graph2):
            return True
    return False


cdef int j_m(ndarray[int, ndim=1] multiple_link_num):
    """Return value of Jm."""
    cdef int num
    cdef int i = 3
    cdef double c = 0
    for num in multiple_link_num[1:]:
        c += i / 2 * num
        i += 1
    return <int>c


cdef int jp_m(ndarray[int, ndim=1] multiple_link_num):
    """Return value of J'm."""
    # Number of multiple links.
    cdef int n_m = sum(multiple_link_num[1:])
    if n_m <= 1:
        return 0
    elif n_m % 2 == 0:
        return <int>((3 * (n_m - 1) - 1) / 2)
    else:
        return <int>((3 * (n_m - 1) - 2) / 2)


cdef tuple n_c(ndarray[int, ndim=1] multiple_link_num):
    """Return all values of Nc."""
    cdef int j_m_v = j_m(multiple_link_num)
    cdef int jp_m_v = jp_m(multiple_link_num)
    return max(1, j_m_v - jp_m_v), min(multiple_link_num[0], j_m_v)


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

    # Number of all joints.
    cdef int joint_count = sum(link_num_)

    # NumPy array type.
    cdef ndarray[int, ndim=1] link_num = np_array(link_num_)

    # Contracted link.
    cdef int n_c_max, n_c_min
    n_c_max, n_c_min = n_c(link_num)
    if n_c_max < n_c_min:
        n_c_max, n_c_min = n_c_min, n_c_max
    # NL2 - Nc + 2
    cdef int i_max = link_num[0] - n_c_min + 2
    cdef int i_min = link_num[0] - n_c_max + 2
    print(i_max, i_min)

    """
    # Number of multiple link joints.
    cdef int mj = 0
    cdef int i
    for i in range(1, len(link_num_)):
        mj += (i + 2) * link_num_[i]

    # Contracted link.
    cdef int cj = link_num_[0]
    cdef int s
    cdef tuple m
    for m in product(range(cj + 1), repeat=cj):
        s = 0
        for i in range(cj):
            s += (i + 1) * m[i]
        if s == cj:
            print(m)
    """

    cdef list answer = []
    # Return graph list and time.
    return answer, (time() - t0)
