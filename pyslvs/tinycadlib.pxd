# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from .triangulation cimport EStack
from .expression cimport Coordinate

cdef double radians(double degree) nogil
cpdef Coordinate plap(Coordinate c1, double d0, double a0, Coordinate c2=*, bint inverse=*)
cpdef Coordinate pllp(Coordinate c1, double d0, double d1, Coordinate c2, bint inverse=*)
cpdef Coordinate plpp(Coordinate c1, double d0, Coordinate c2, Coordinate c3, bint inverse=*)
cpdef Coordinate pxy(Coordinate c1, double x, double y)

cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints, object angles=*)
