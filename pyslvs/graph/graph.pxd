# -*- coding: utf-8 -*-
# cython: language_level=3

"""Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.map cimport map

ctypedef unsigned long long ullong
ctypedef map[int, int] imap

cpdef list link_assortment(Graph g)
cpdef list contracted_link_assortment(Graph g)

cdef class Graph:
    # Graph(edges)

    cdef readonly tuple edges
    cdef readonly tuple vertices
    cdef dict adj

    cpdef void add_vertices(self, object new_nodes)
    cpdef void add_edge(self, int n1, int n2)
    cpdef void add_path(self, object new_nodes)
    cpdef void remove_edge(self, int n1, int n2)

    cpdef int dof(self)
    cpdef tuple neighbors(self, int n)
    cpdef dict degrees(self)
    cpdef ullong degree_code(self)
    cpdef double[:, :] adjacency_matrix(self)

    cpdef bint is_connected(self, int without=*)
    cpdef bint has_cut_link(self)
    cpdef bint is_degenerate(self)
    cpdef bint has_triangle(self)
    cpdef bint is_isomorphic(self, Graph g)
    cpdef bint is_isomorphic_vf2(self, Graph g)
    cpdef bint is_isomorphic_degree_code(self, Graph g)

    cpdef Graph duplicate(self, object vertices, int times)
    cpdef Graph copy(self)
