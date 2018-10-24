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

    cdef readonly tuple edges
    cdef tuple nodes
    cdef dict adj

    cdef inline tuple neighbors(self, int n)
    cpdef int dof(self)
    cpdef bool is_connected(self, int with_out=*)
    cpdef bool has_cut_link(self)
    cpdef bool is_degenerate(self)
    cpdef bool is_isomorphic(self, Graph graph)
    cdef list link_types(self)
    cdef int number_of_edges(self, int u, int v)
    cdef list multi_contracted_links(self)
    cdef bool has_triangles(self)
    cpdef Graph copy(self)
