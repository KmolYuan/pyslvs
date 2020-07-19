# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from .topo_config cimport EStack
from .expression cimport Coord

cdef double radians(double degree) nogil

cpdef Coord pxy(Coord c1, double x, double y)
cpdef Coord ppp(Coord c1, Coord c2, Coord c3)
cpdef Coord plap(Coord c1, double d0, double a0, Coord c2=*, bint inverse=*)
cpdef Coord pllp(Coord c1, double d0, double d1, Coord c2, bint inverse=*)
cpdef Coord plpp(Coord c1, double d0, Coord c2, Coord c3, bint inverse=*)
cpdef Coord palp(Coord c1, double a0, double d0, Coord c2, bint inverse=*)

cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints,
                        object angles=*)
