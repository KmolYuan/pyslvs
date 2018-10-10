# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from cpython cimport bool
from numpy cimport ndarray


cdef class Graph:

    """NetworkX-like graph class."""

    cdef public tuple edges
    cdef tuple nodes
    cdef dict adj

    cdef inline tuple neighbors(self, int n)
    cdef Graph compose(self, Graph graph)
    cdef bool out_of_limit(self, ndarray limit)
    cpdef bool has_triangles(self)
    cpdef bool is_connected(self)
    cpdef bool is_isomorphic(self, Graph graph)
    cdef list links(self)
    cdef int number_of_edges(self, int u, int v)
