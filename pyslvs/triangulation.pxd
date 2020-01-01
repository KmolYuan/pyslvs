# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.list cimport list as clist
from libcpp.pair cimport pair

ctypedef pair[int, int] symbol

cdef str symbol_str(symbol p)
cpdef ExpressionStack vpoints_configure(object vpoints_, object inputs, dict status=*)

cdef enum Label:
    P_LABEL
    L_LABEL
    A_LABEL
    S_LABEL

cdef enum Func:
    PLA
    PLAP
    PLLP
    PLPP
    PXY

cdef struct Expression:
    bint op
    Func func
    symbol v1
    symbol v2
    symbol c1
    symbol c2
    symbol c3
    symbol c4

cdef class ExpressionStack:

    cdef clist[Expression] stack

    cdef void add_pla(self, symbol c1, symbol v1, symbol v2, symbol target)
    cdef void add_plap(self, symbol c1, symbol v1, symbol v2, symbol c2, symbol target)
    cdef void add_pllp(self, symbol c1, symbol v1, symbol v2, symbol c2, symbol target)
    cdef void add_plpp(self, symbol c1, symbol v1, symbol c2, symbol c3, symbol target, bint op)
    cdef void add_pxy(self, symbol c1, symbol v1, symbol v2, symbol target)

    cpdef list as_list(self)
