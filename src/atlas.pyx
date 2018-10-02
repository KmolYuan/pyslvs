# -*- coding: utf-8 -*-
# cython: language_level=3

"""Structure synthesis."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

cimport cython
from typing import (
    Tuple,
    Iterator,
)
from itertools import combinations, product
from time import time
from cpython cimport bool
from numpy cimport ndarray
from numpy import (
    zeros as np_zeros,
    int32 as np_int32,
)
from graph cimport Graph


cdef inline bool is_isomorphic(Graph graph1, list answer):
    """Return True if the graph is isomorphic with list."""
    cdef Graph graph2
    for i, graph2 in enumerate(answer):
        if graph1.is_isomorphic(graph2):
            return True
    return False


cdef inline list get_bind(int i, tuple connection):
    cdef tuple c
    return [c for c in connection if (i in c)]


def _graphs(c: Iterator[Tuple[Tuple[int, int], ...]]) -> Iterator[Graph]:
    """Graph generator."""
    cdef tuple m
    for m in c:
        yield Graph(m)


cdef inline int factorial(int x):
    """Factorial function."""
    cdef int m = x
    cdef int i
    if x <= 1:
        return 1
    else:
        for i in range(1, x):
            m *= i
        return m


@cython.cdivision(True)
cdef int c_count(int n, int k):
    """Combination count."""
    return factorial(n) / (factorial(k) * factorial(n - k))


cpdef tuple topo(
    object link_num,
    bool degenerate = True,
    object job_func = None,
    object stop_func = None
):
    """Linkage mechanism topological function.
    
    link_num = [L2, L3, L4, ...]
    links = [[number_code]: joint_number, ...]
    """
    # Initial time.
    cdef double t0 = time()
    
    # Number of all joints.
    cdef int joint_count = sum(link_num)
    
    # Number of multiple link joints.
    cdef int mj = 0
    cdef int i
    for i in range(1, len(link_num)):
        mj += (i + 2) * link_num[i]
    
    # Contracted link.
    cdef int cj = link_num[0]
    cdef int s
    cdef tuple m
    for m in product(range(cj + 1), repeat=cj):
        s = 0
        for i in range(cj):
            s += (i + 1) * m[i]
        if sum(m) * 2 >= mj:
            continue
        if s == cj:
            print(m)
    
    # Number of joins in each links.
    cdef ndarray links = np_zeros((joint_count,), dtype=np_int32)
    cdef int j, t, name
    for i in range(joint_count):
        name = i
        for j, t in enumerate(link_num):
            if i < t:
                links[name] = j + 2
                break
            i -= t
        else:
            raise RuntimeError("Invalid number")
    
    # binds = [(0, 1), (0, 2), ..., (1, 2), (1, 3), ...]
    cdef tuple binds = tuple(combinations(range(joint_count), 2))
    
    # Result list.
    cdef list edges_combinations = []
    # Iterator.
    cdef object match, matched
    cdef list current_binds
    cdef int match_count
    
    # Find ALL results.
    cdef int link, count, progress_value
    cdef Graph graph1, graph2, graph3
    for link, count in enumerate(links):
        # Other of joints that the link connect with.
        current_binds = get_bind(link, binds)
        match = _graphs(combinations(current_binds, count))
        match_count = c_count(len(current_binds), count)
        
        # First population.
        if not edges_combinations:
            edges_combinations.extend(match)
            continue
        
        if job_func:
            progress_value = len(edges_combinations) * match_count
            job_func(
                f"Match link # {link} / {len(links) - 1}\n"
                f"Connections: {count}\n\n"
                f"Feasible: {len(edges_combinations)}\n"
                f"Matching: {match_count}\n"
                f"Possibilities: {progress_value}",
                progress_value
            )
        
        # Collecting.
        matched = product(edges_combinations, match)
        edges_combinations.clear()
        for graph1, graph2 in matched:
            if stop_func and stop_func():
                break
            graph3 = graph1.compose(graph2)
            
            # Out of limit.
            if graph3.out_of_limit(links):
                continue
            
            # Has triangles.
            if degenerate and graph3.has_triangles():
                continue
            
            edges_combinations.append(graph3)
    
    if job_func:
        progress_value = len(edges_combinations)
        job_func(
            f"Verify the graphs ...\n"
            f"Count: {progress_value}",
            progress_value
        )
    
    cdef list answer = []
    for graph1 in edges_combinations:
        if stop_func and stop_func():
            break
        if not graph1.is_connected():
            continue
        if is_isomorphic(graph1, answer):
            continue
        answer.append(graph1)
    
    # Return graph list and time.
    return answer, (time() - t0)
