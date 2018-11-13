# -*- coding: utf-8 -*-
# cython: language_level=3

"""Number synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from numpy import (
    int16,
    array as np_array,
    zeros as np_zeros,
)
from numpy cimport (
    ndarray,
    int16_t,
)


cdef inline list product(object pool, int repeat):
    """Product function as same as iteration tools."""
    cdef int i, y
    cdef list x, tmp_list
    cdef list result = [[]]
    for i in range(repeat):
        tmp_list = []
        for x in result:
            for y in pool:
                tmp_list.append(x + [y])
        result = tmp_list
    return result


cdef inline int m_max(int nl, int nj):
    """Find max number of joint on each link.
    
    + nl <= nj and nj <= (2 * nl - 3)
    + (2 * nl - 3) <= nj and nj <= (nl * (nl - 1) / 2)
    + other exceptions (return -1).
    """
    if nl <= nj <= (2 * nl - 3):
        return nj - nl + 2
    if nl == nj == 0:
        return -1
    if (2 * nl - 3) <= nj <= (nl * (nl - 1) / 2):
        return nl - 1
    return -1


cdef inline int sum_factors(list factors):
    """F0*N2 + F1*N3 + F2*N4 + ... + Fn*N(n+2)"""
    cdef int factor = 0
    cdef int i
    for i in range(len(factors)):
        factor += factors[i] * (i + 2)
    return factor


cpdef tuple number_synthesis(int nl, int nj):
    """Number synthesis try-error function."""
    cdef list result = []
    cdef int m_max_v = m_max(nl, nj)
    if m_max_v == -1:
        raise Exception("incorrect mechanism.")
    cdef int i, p
    cdef list symbols
    for symbols in product(range(nl + 1), m_max_v - 2):
        nl_m_max = nl - sum(symbols)
        if nl_m_max < 0:
            continue
        symbols.append(nl_m_max)
        if sum_factors(symbols) == (2 * nj):
            result.append(symbols)
    return tuple(result)


cdef int j_m(int16_t[:] link_num):
    """Return value of Jm."""
    cdef int num
    cdef int i = 3
    cdef float c = 0
    for num in link_num[1:]:
        c += i / 2 * num
        i += 1
    return <int>c


cdef int j_m_p(int16_t[:] link_num):
    """Return value of Jm'. This is improved function.
    
    + Origin equation:
    if n_m % 2 == 0:
        return <int>((3 * (n_m - 1) - 1) / 2)
    else:
        return <int>((3 * (n_m - 1) - 2) / 2)
    """
    # Number of multiple links.
    cdef int n_m = sum(link_num[1:])
    if n_m <= 1:
        return 0
    elif n_m == 2:
        return 1
    else:
        return 3 * (n_m - 2)


cdef tuple n_c(int16_t[:] link_num):
    """Return all values of Nc."""
    cdef int j_m_v = j_m(link_num)
    cdef int j_m_p_v = j_m_p(link_num)
    return max(1, j_m_v - j_m_p_v), min(link_num[0], j_m_v)


cpdef list contracted_link(list link_num_list):
    """Generate the contracted link assortments."""
    cdef ndarray[int16_t, ndim=1] link_num

    if len(link_num_list) == 1:
        link_num = np_zeros(link_num_list[0], dtype=int16)
        link_num[-1] = 1
        return [tuple(link_num.tolist())]

    link_num = np_array(link_num_list, ndmin=1, dtype=int16)

    # Contracted link.
    cdef int n_c_min, n_c_max
    n_c_min, n_c_max = n_c(link_num)

    # NL2 - Nc + 2
    cdef int i_max = link_num[0] - n_c_min + 2

    # Matching formula.
    cdef int count, factor, index
    cdef list m
    cdef list cj_list = []
    for m in product(range(link_num[0] + 1), repeat=i_max):
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
            cj_list.append(tuple(m))

    return cj_list
