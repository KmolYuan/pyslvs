# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from .tinycadlib cimport Expr

ctypedef pair[int, int] Sym

cpdef EStack t_config(object vpoints_, object inputs_, object status_=*)

cdef class EStack:
    cdef vector[Expr] stack

    cdef void add_pxy(self, Sym c1, Sym v1, Sym v2, Sym t) nogil
    cdef void add_ppp(self, Sym c1, Sym c2, Sym c3, Sym t) nogil
    cdef void add_pla(self, Sym c1, Sym v1, Sym v2, Sym t) nogil
    cdef void add_plap(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t) nogil
    cdef void add_pllp(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t) nogil
    cdef void add_plpp(self, Sym c1, Sym v1, Sym c2, Sym c3, Sym t, bint op) nogil
    cdef void add_palp(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t, bint op) nogil
    cpdef list as_list(self)
