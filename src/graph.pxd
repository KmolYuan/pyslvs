# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.map cimport map as c_map

ctypedef c_map[int, int] c_map_int


cdef class Graph:

    """NetworkX-like graph class."""

    cdef readonly tuple edges
    cdef readonly tuple nodes
    cdef dict adj

    cpdef void add_edge(self, int n1, int n2)
    cpdef void add_nodes_from(self, tuple nodes)

    cpdef int dof(self)
    cpdef tuple neighbors(self, int n)
    cdef c_map_int degrees(self)

    cpdef bint is_connected(self, int with_out=*)
    cpdef bint has_cut_link(self)
    cpdef bint is_degenerate(self)
    cpdef bint is_isomorphic(self, Graph graph)
    cdef bint has_triangles(self)

    cpdef Graph copy(self)
