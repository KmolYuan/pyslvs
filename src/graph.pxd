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

    cpdef void add_edge(self, int n1, int n2)
    cpdef void add_nodes_from(self, tuple nodes)
    cdef tuple neighbors(self, int n)
    cpdef list degrees(self)
    cpdef int dof(self)
    cpdef bool is_connected(self, int with_out=*)
    cpdef bool has_cut_link(self)
    cpdef bool is_degenerate(self)
    cpdef bool is_isomorphic(self, Graph graph)
    cpdef bool is_planar(self)
    cdef list link_types(self)
    cdef int node_distance(self, int u, int v)
    cdef list multi_contracted_links(self)
    cdef bool has_triangles(self)
    cpdef Graph copy(self)
    cpdef Graph subgraph(self, tuple nodes)
