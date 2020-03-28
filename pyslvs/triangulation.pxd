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

ctypedef pair[int, int] sym

cdef str symbol_str(sym p)
cpdef EStack t_config(object vpoints_, object inputs, object status=*)

cdef enum Label:
    P_LABEL
    L_LABEL
    I_LABEL
    S_LABEL

cdef enum Func:
    PLA
    PLAP
    PLLP
    PLPP
    PXY

cdef struct Expr:
    bint op
    Func func
    sym v1
    sym v2
    sym c1
    sym c2
    sym c3
    sym c4

cdef class EStack:

    cdef vector[Expr] stack

    cdef void add_pla(self, sym c1, sym v1, sym v2, sym target) nogil
    cdef void add_plap(self, sym c1, sym v1, sym v2, sym c2, sym target) nogil
    cdef void add_pllp(self, sym c1, sym v1, sym v2, sym c2, sym target) nogil
    cdef void add_plpp(self, sym c1, sym v1, sym c2, sym c3, sym target, bint op) nogil
    cdef void add_pxy(self, sym c1, sym v1, sym v2, sym target) nogil
    cpdef list as_list(self)
