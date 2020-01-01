# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True, cdivision=True

"""Number synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import (
    int16 as np_int,
    array as np_array,
    zeros as np_zeros,
    sum as np_sum,
    prod as np_prod,
    repeat as np_repeat,
    arange,
)


@cython.boundscheck(False)
cdef int16_t[:, :] product(tuple pool, object stop_func):
    """Product function as same as iteration tools.

    The pool is created by range(n).
    """
    if not pool:
        return np_array([], dtype=np_int)

    cdef int16_t[:] tmp1
    cdef int16_t[:, :] tmp2
    cdef int16_t[:] array0 = arange(pool[0], dtype=np_int)

    cdef int n = np_prod(pool)
    cdef int16_t[:, :] out = np_zeros((n, len(pool)), dtype=np_int)

    cdef int array0_size = len(array0)
    cdef int m = n / array0_size
    tmp1 = np_repeat(array0, m)
    out[:, 0] = tmp1

    cdef int j
    if pool[1:]:
        tmp2 = product(pool[1:], stop_func)
        out[0:m, 1:] = tmp2
        for j in range(1, array0_size):
            if stop_func is not None and stop_func():
                return out
            out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
    return out


cdef int16_t[:, :] _product(int pool_size, int repeat, object stop_func):
    return product((pool_size,) * repeat + (1,), stop_func)


cdef inline int _m_max(int nl, int nj) nogil:
    """Find max number of joint on each link.

    + nl <= nj and nj <= (2 * nl - 3)
    + (2 * nl - 3) <= nj and nj <= (nl * (nl - 1) / 2)
    + other exceptions (return -1).
    """
    if nl <= nj <= (2 * nl - 3):
        return nj - nl + 2
    if nl == nj == 0:
        return -1
    if 2 * nl - 3 <= nj <= <double>(nl * (nl - 1)) / 2:
        return nl - 1
    return -1


cdef inline int _sum_factors(int16_t[:] factors):
    """F0*N2 + F1*N3 + F2*N4 + ... + Fn*N(n+2)"""
    cdef int factor = 0
    cdef int i
    for i in range(len(factors)):
        factor += factors[i] * (i + 2)
    return factor


cpdef list link_synthesis(int nl, int nj, object stop_func = None):
    """Number synthesis try and error function."""
    result = []
    cdef int m_max_v = _m_max(nl, nj)
    if m_max_v == -1:
        raise ValueError("incorrect mechanism.")

    cdef int i, p
    cdef int16_t[:] symbols
    for symbols in _product(nl + 1, m_max_v - 2, stop_func):
        nl_m_max = nl - np_sum(symbols)
        if nl_m_max < 0:
            continue
        symbols[-1] = nl_m_max
        if _sum_factors(symbols) == nj * 2:
            result.append(tuple(symbols))
    return result


cdef inline int _j_m(int16_t[:] link_num):
    """Return value of JM."""
    cdef int num
    cdef int i = 3
    cdef float c = 0
    for num in link_num[1:]:
        c += <double>i / 2 * num
        i += 1
    return <int>c


cdef inline int _j_m_p(int n_m) nogil:
    """Return value of J'M. This is improved function.

    + Origin equation:
    if n_m % 2 == 0:
        return <int>((3 * (n_m - 1) - 1) / 2)
    else:
        return <int>((3 * (n_m - 1) - 2) / 2)
    """
    # Number of multiple links.
    if n_m <= 1:
        return 0
    elif n_m == 2:
        return 1
    else:
        return 3 * (n_m - 2)


cpdef list contracted_link_synthesis(object link_num_list, object stop_func = None):
    """Generate the contracted link assortment."""
    cdef int16_t[:] link_num
    if len(link_num_list) == 1:
        link_num = np_zeros(link_num_list[0], dtype=np_int)
        link_num[-1] = 1
        return [tuple(link_num)]

    link_num = np_array(link_num_list, dtype=np_int)

    # Contracted link
    cdef int j_m_v = _j_m(link_num)
    cdef int n_c_min = max(1, j_m_v - _j_m_p(sum(link_num[1:])))
    cdef int n_c_max = min(link_num[0], j_m_v)

    # i = NL2 - NC + 2
    cdef int i_max = min(link_num[0], link_num[0] - n_c_min + 2)

    # Matching formula
    cdef int count, factor, index
    cdef float last_factor
    cdef int16_t[:] m
    cj_list = []
    for m in _product(link_num[0] + 1, i_max - 1, stop_func):
        count = 0
        index = 1
        for factor in m[:-1]:
            count +=  factor * index
            index += 1

        # Check if the last factor is a natural number.
        last_factor = <double>(link_num[0] - count) / index
        factor = <int>last_factor
        if last_factor < 0 or last_factor != factor:
            continue

        m[-1] = factor
        if n_c_min <= np_sum(m) <= n_c_max:
            cj_list.append(tuple(m))

    return cj_list
