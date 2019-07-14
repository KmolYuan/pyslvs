# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.map cimport map as cmap

ctypedef cmap[int, int] imap


cpdef list link_assortment(Graph g)
cpdef list contracted_link_assortment(Graph g)


cdef class Graph:

    # Graph(edges)

    cdef readonly tuple edges
    cdef readonly tuple nodes
    cdef dict adj

    cpdef void add_nodes(self, object new_nodes)
    cpdef void add_edge(self, int n1, int n2)
    cpdef void add_path(self, object new_nodes)
    cpdef void remove_edge(self, int n1, int n2)

    cpdef int dof(self)
    cpdef tuple neighbors(self, int n)
    cdef dict degrees(self)

    cpdef bint is_connected(self, int without=*)
    cpdef bint has_cut_link(self)
    cpdef bint is_degenerate(self)
    cpdef bint is_isomorphic(self, Graph graph)
    cdef bint has_triangles(self)

    cpdef Graph duplicate(self, object nodes, int times)
    cpdef Graph copy(self)
