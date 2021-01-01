# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structural synthesis.

The algorithm reference:
+ On the Number Synthesis of Kinematic Chains
    + author: Hong-Sen Yan, Yu-Ting Chiu
    + Mechanism and Machine Theory, Volume 89, 2015, Pages 128-144, ISSN 0094-114X
    + https://doi.org/10.1016/j.mechmachtheory.2014.08.012

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from time import process_time
from logging import getLogger
from collections import Counter
cimport cython
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.map cimport map
from numpy import (
    int16 as i16,
    array as np_array,
    zeros as np_zeros,
    ones as np_ones,
    sum as np_sum,
    subtract as np_sub,
    multiply as np_mul,
    floor_divide as np_div,
    any as np_any,
    prod as np_prod,
    repeat as np_repeat,
    arange,
)
from .graph cimport Graph, link_assortment
from .planar cimport is_planar

ctypedef unsigned int uint
ctypedef unsigned long long ullong
ctypedef pair[int, int] ipair
ctypedef map[int, int] imap

logger = getLogger()


@cython.boundscheck(False)
@cython.wraparound(False)
cdef short[:, :] _product(tuple pool, object stop_func):
    """Product function as same as iteration tools.

    The pool is created by range(n).
    """
    if not pool:
        return np_array([], dtype=i16)
    cdef short[:] tmp1
    cdef short[:, :] tmp2
    cdef short[:] array0 = arange(pool[0], dtype=i16)
    cdef int n = np_prod(pool)
    cdef short[:, :] out = np_zeros((n, len(pool)), dtype=i16)
    cdef int array0_size = len(array0)
    cdef int m = n / array0_size
    tmp1 = np_repeat(array0, m)
    out[:, 0] = tmp1
    cdef int j
    if pool[1:]:
        tmp2 = _product(pool[1:], stop_func)
        out[0:m, 1:] = tmp2
        for j in range(1, array0_size):
            if stop_func is not None and stop_func():
                return out
            out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
    return out


cdef short[:, :] _simple_product(int pool_size, int repeat, object stop_func):
    return _product((pool_size,) * repeat + (1,), stop_func)


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


cdef inline int _sum_factors(short[:] factors):
    """F0*N2 + F1*N3 + F2*N4 + ... + Fn*N(n+2)"""
    cdef int factor = 0
    cdef int i
    for i in range(len(factors)):
        factor += factors[i] * (i + 2)
    return factor


cpdef list link_synthesis(int nl, int nj, object stop_func = None):
    """Return link assortment by number of links `nl` and number of joints `nj`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    result = []
    cdef int m_max_v = _m_max(nl, nj)
    if m_max_v == -1:
        raise ValueError("incorrect mechanism.")

    cdef int i, p
    cdef short[:] symbols
    for symbols in _simple_product(nl + 1, m_max_v - 2, stop_func):
        nl_m_max = nl - np_sum(symbols)
        if nl_m_max < 0:
            continue
        symbols[-1] = nl_m_max
        if _sum_factors(symbols) == nj * 2:
            result.append(tuple(symbols))
    return result


cdef inline int _j_m(short[:] link_num):
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
    """Return contracted link assortment by link assortment `link_num_list`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    cdef short[:] link_num
    if len(link_num_list) == 1:
        link_num = np_zeros(link_num_list[0], dtype=i16)
        link_num[-1] = 1
        return [tuple(link_num)]

    link_num = np_array(link_num_list, dtype=i16)

    # Contracted link
    cdef int j_m_v = _j_m(link_num)
    cdef int n_c_min = max(1, j_m_v - _j_m_p(sum(link_num[1:])))
    cdef int n_c_max = min(link_num[0], j_m_v)

    # i = NL2 - NC + 2
    cdef int i_max = min(link_num[0], link_num[0] - n_c_min + 2)

    # Matching formula
    cdef int count, factor, index
    cdef float last_factor
    cdef short[:] m
    cj_list = []
    for m in _simple_product(link_num[0] + 1, i_max - 1, stop_func):
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


cdef inline short[:] _nonzero_index(short[:] array):
    """Return the number of nonzero numbers."""
    counter = []
    cdef int i, n
    for i, n in enumerate(array):
        if n != 0:
            counter.append(i)
    return np_array(counter, dtype=i16)


cdef inline ullong _factorial(int n):
    """Return (n!)."""
    cdef ullong ret = 1
    cdef int i
    for i in range(2, n + 1):
        ret *= i
    return ret


cdef inline int _gcd(int a, int b) nogil:
    """Return Greatest Common Divisor of a and b.

    Only for positive numbers.
    """
    cdef int r
    while b > 0:
        r = a % b
        a = b
        b = r
    return a


cdef inline short[:] _labels(short[:] numbers, int index, int offset):
    """Generate labels from numbers."""
    cdef int i, num
    labels = []
    for num in numbers[offset:]:
        for i in range(num):
            labels.append(index)
        index += 1
    return np_array(labels, dtype=i16)


cdef inline Graph _multigraph(short[:] counter, int n):
    """Get multigraph from n x n matrix."""
    edges = {}
    cdef int c = 0
    cdef int i, j
    for i in range(n):
        for j in range(n):
            if i >= j:
                continue
            edges[i, j] = counter[c]
            c += 1
    return Graph(Counter(edges).elements())


@cython.boundscheck(False)
cdef inline void _gauss_elimination(
    list result,
    short[:] limit,
    short[:, :] f_matrix,
    int n,
    int var_count
):
    """Gauss elimination for (n, n + 1) matrix."""
    cdef int i, j, d
    cdef short[:] tmp1, tmp2
    for j in range(var_count):
        # Remove all coefficients of index [i] to zero
        for i in range(n):
            if f_matrix[i, j] != 0 and not np_any(f_matrix[i, :j]):
                d = i
                tmp2 = np_div(f_matrix[i, :], f_matrix[i, j])
                break
        else:
            continue

        for i in range(n):
            if i == d or f_matrix[i, j] == 0:
                continue
            tmp1 = np_sub(f_matrix[i, :], np_mul(tmp2, f_matrix[i, j]))
            f_matrix[i, :] = tmp1

    # Answer
    cdef short[:] answer = -np_ones(var_count, dtype=i16)

    # Determined solution
    cdef int c, k
    for i in range(n):
        c = 0
        for j in range(var_count):
            # Derivation (has answer)
            if answer[j] != -1 and f_matrix[i, j] != 0:
                f_matrix[i, -1] -= f_matrix[i, j] * answer[j]
                f_matrix[i, j] = 0

            # Nonzero coefficient
            if f_matrix[i, j] != 0:
                d = j
                c += 1

        if c != 1:
            continue

        j = d
        k = f_matrix[i, j]
        c = f_matrix[i, -1]
        if k != 1:
            if k < 0:
                k = -k
                c = -c
            if c < 0:
                return
            d = _gcd(k, c)
            k /= d
            c /= d
        if c < 0:
            return
        answer[j] = c

    # Result
    _test_contracted_graph(_multigraph(answer, n), result)


@cython.boundscheck(False)
cdef void _nest_do(
    list result,
    short[:] answer,
    short[:, :] f_matrix,
    int i,
    int n,
    object stop_func
):
    """Nest do loop."""
    if i >= n:
        # Result
        if <int>np_sum(answer) == <int>np_sum(f_matrix[:, -1]) / 2:
            _test_contracted_graph(_multigraph(answer, n), result)
        return

    cdef short[:] coefficients = _nonzero_index(f_matrix[i, :-1])
    cdef int c = 0
    cdef int j = -1
    cdef int d
    for d in coefficients:
        if answer[d] == -1:
            c += 1
            j = d

    if c == 0:
        if <int>np_sum(np_mul(answer, f_matrix[i, :-1])) == f_matrix[i, -1]:
            _nest_do(result, answer, f_matrix, i + 1, n, stop_func)
        return
    elif c == 1:
        c = 0
        for d in coefficients:
            if answer[d] != -1:
                c += answer[d] * f_matrix[i, d]
        answer[j] = (f_matrix[i, -1] - c) / f_matrix[i, j]
        _nest_do(result, answer, f_matrix, i + 1, n, stop_func)
        return

    cdef short[:] combine, answer_copy
    for combine in _product((f_matrix[i, -1],) * c, stop_func):
        if stop_func is not None and stop_func():
            return

        c = 0
        d = 0
        for j in coefficients:
            if answer[j] == -1:
                c += combine[d] * f_matrix[i, j]
                d += 1
            else:
                c += answer[j] * f_matrix[i, j]
        if c != f_matrix[i, -1]:
            continue

        answer_copy = answer.copy()
        d = 0
        for j in coefficients:
            # Pass to answer
            if answer_copy[j] == -1:
                answer_copy[j] = combine[d]
                d += 1
        _nest_do(result, answer_copy, f_matrix, i + 1, n, stop_func)


cdef inline bint _is_isomorphic(Graph g, list result):
    """Return true if graph is isomorphic with result list."""
    cdef Graph h
    for h in result:
        if g.is_isomorphic(h):
            return True
    return False


cdef inline void _test_contracted_graph(Graph g, list result):
    """Test the contracted graph."""
    if not g.edges:
        return
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


@cython.boundscheck(False)
cdef inline void _contracted_graph(
    list result,
    short[:] limit,
    object stop_func
):
    """Synthesis of contracted graphs."""
    cdef int n = len(limit)
    cdef int var_count = n * (n - 1) / 2
    cdef short[:, :] f_matrix = np_zeros((n, var_count + 1), dtype=i16)
    f_matrix[:, -1] = limit

    # Equations
    cdef int i, j, k, c
    for i in range(n):
        c = 0
        for j in range(n):
            for k in range(n):
                if j >= k:
                    continue
                if i in {j, k}:
                    f_matrix[i, c] = 1
                c += 1

    if n >= var_count:
        # Fast solution by Gauss elimination.
        _gauss_elimination(result, limit, f_matrix, n, var_count)
    else:
        # Nest do loop method.
        _nest_do(result, -np_ones(var_count, dtype=i16), f_matrix, 0, n, stop_func)


cdef inline void _dyad_insert(Graph g, frozenset edge, int amount):
    """Insert dyad to the graph."""
    if amount < 1:
        return

    cdef int last_num = max(g.vertices) + 1
    cdef int n1, n2
    n1, n2 = edge

    path = [n1]
    path.extend(range(last_num, last_num + amount))
    path.append(n2)

    g.remove_edge(n1, n2)
    g.add_path(path)


cdef inline void _permute_combine(
    short[:] limit,
    list combine_list,
    tuple pick_list,
    object stop_func
):
    """Permutation of combined list."""
    cdef int n = len(limit)
    if n < 1:
        return

    cdef vector[int] indices = range(n)
    cdef short[:] cycles = np_zeros(n, dtype=i16)
    cdef short[:] pool = limit
    permute_list = set()
    cdef int i, j
    for i, j in enumerate(range(n, 0, -1)):
        cycles[i] = j
    permute_list.add(tuple(pool[indices[i]] for i in range(n)))
    cdef vector[int].iterator it2 = indices.begin()

    cdef int tmp
    while True:
        if stop_func is not None and stop_func():
            break

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
    for tmp_array in permute_list:
        combine_list.append(tuple(zip(pick_list, tmp_array)))


cdef inline list _contracted_links(tuple edges, short[:] limit, object stop_func):
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
    pool_list = tuple([frozenset(edge) for edge in edges])
    pool = Counter(pool_list)
    # The list including required edge(s).
    confirm_list = pool - Counter(set(pool_list))
    # Simplified the pool
    pool_list = tuple((pool - confirm_list).elements())
    cdef int confirm_size = sum(confirm_list.values())
    if confirm_size > pick_count or pick_count > len(edges):
        return []

    pick_count -= confirm_size
    cdef int pool_size = len(pool_list)
    cdef short[:] indices = np_zeros(pick_count, dtype=i16)
    cdef int i
    for i in range(pick_count):
        indices[i] = i
    pick_list = Counter()
    combine_set = set()
    # Combinations loop with number checking.
    cdef int n1, n2
    while True:
        if stop_func is not None and stop_func():
            break

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
    combine_list = []
    pool_list = tuple(confirm_list.elements())
    for picked in combine_set:
        _permute_combine(limit, combine_list, pool_list + picked, stop_func)
    return combine_list


cdef inline void _graph_atlas(
    list result,
    list contracted_graph,
    short[:] limit,
    uint no_degenerate,
    object stop_func
):
    """Synthesis of atlas."""
    cdef int n
    cdef Graph cg, g
    for cg in contracted_graph:
        for combine in _contracted_links(cg.edges, limit, stop_func):
            g = Graph.__new__(Graph, cg.edges)
            for edge, n in combine:
                _dyad_insert(g, edge, n)
            _test_graph(g, result, no_degenerate)


cdef inline list _loop_chain(int num):
    """Loop chain of num."""
    cdef int i
    cdef int b = 0
    chain = []
    for i in range(1, num):
        chain.append((b, i))
        b = i
    chain.append((0, b))
    return chain


cpdef list contracted_graph(object link_num_list, object stop_func = None):
    """Generate contracted graphs by link assortment `link_num`.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    if not link_num_list:
        return []

    # Initial time
    cdef double t0 = process_time()
    cdef short[:] link_num = np_array(link_num_list, ndmin=1, dtype=i16)
    logger.debug(f"Link assortment: {list(link_num)}")

    # Multiple links
    cdef short[:] m_link = _labels(link_num, 3, 1)

    cdef imap m_limit, count
    cdef int num
    cdef int i = 0
    for num in m_link:
        m_limit[i] = num
        count[i] = 0
        i += 1

    # Synthesis of contracted graphs
    cg_list = []
    _contracted_graph(cg_list, m_link, stop_func)
    logger.debug(f"Contracted graph(s): {len(cg_list)}, time: {process_time() - t0}")
    return cg_list


cpdef list conventional_graph(
    list cg_list,
    object c_j_list,
    uint no_degenerate = 1,
    object stop_func = None
):
    """Generate conventional graphs by contracted graphs `cg_list` and
    contracted link assortment `c_j_list`.

    The degenerate setting `no_degenerate` has following option:

    + `0`: No degenerate.
    + `1`: Only degenerate.
    + Else: All graphs.

    The check stop function `stop_func` object for GUI or subprocess,
    return `True` to terminate this function.
    """
    # Initial time
    cdef double t0 = process_time()
    logger.debug(f"Contracted link assortment: {list(c_j_list)}")

    # Synthesis of contracted link and multiple link combination.
    cdef short[:] c_j = np_array(c_j_list, ndmin=1, dtype=i16)
    result = []
    cdef int i, num
    if not cg_list:
        if 1 not in c_j_list:
            return []
        # Single loop - ring graph (special case)
        i = 1
        for num in c_j:
            if num == 1:
                break
            i += 1
        else:
            raise ValueError("Invalid assortment!")
        result.append(Graph.__new__(Graph, _loop_chain(i)))
        logger.debug(f"Count: {len(result)}")
        return result

    # Multiple links
    cdef short[:] m_link = np_array(link_assortment(cg_list[0]), ndmin=1, dtype=i16)
    m_link = _labels(m_link, 3, 1)

    # Synthesis of multiple links
    _graph_atlas(result, cg_list, _labels(c_j, 1, 0), no_degenerate, stop_func)

    # Return graph list and time
    logger.debug(f"Count: {len(result)}, time: {process_time() - t0}")
    return result
