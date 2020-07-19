# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.vector cimport vector
from libcpp.pair cimport pair

# Label, int
ctypedef pair[int, int] Sym

cdef str symbol_str(Sym p)
cpdef EStack t_config(object vpoints_, object inputs_, object status_=*)

cdef enum Label:
    P_LABEL
    L_LABEL
    I_LABEL
    A_LABEL
    S_LABEL

cdef enum Func:
    PXY
    PPP
    PLA
    PLAP
    PLLP
    PLPP
    PALP

cdef struct Expr:
    bint op
    Func func
    Sym v1, v2, c1, c2, c3, target

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
