# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from cpython cimport bool


cdef class Graph:

    """NetworkX-like graph class."""

    cdef public tuple edges
    cdef tuple nodes
    cdef dict adj

    cdef inline tuple neighbors(self, int n)
    cpdef bool has_triangles(self)
    cpdef bool is_connected(self)
    cpdef bool is_isomorphic(self, Graph graph)
    cdef list links(self)
    cdef int number_of_edges(self, int u, int v)
    cpdef Graph copy(self)
