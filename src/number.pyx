# -*- coding: utf-8 -*-
# cython: language_level=3

"""Number synthesis."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"


cdef list product(object pool, int repeat=1):
    """Product function as same as iteration tools."""
    cdef list x, tmp_list
    cdef int i, y
    cdef list result = [[]]
    for i in range(repeat):
        tmp_list = []
        for x in result:
            for y in pool:
                tmp_list.append(x + [y])
        result = tmp_list
    return result


cdef inline int mmax(int nl, int nj):
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
    cdef list result = []
    cdef int m_max = mmax(nl, nj)
    if m_max == -1:
        raise Exception("incorrect mechanism.")
    cdef int i, p
    cdef list symbols
    for symbols in product(range(nl + 1), m_max - 2):
        nl_m_max = nl - sum(symbols)
        if nl_m_max < 0:
            continue
        symbols.append(nl_m_max)
        if sum_factors(symbols) == (2 * nj):
            result.append(symbols)
    return tuple(result)
