# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from cpython cimport bool


cdef class Coordinate:
    cdef readonly double x, y
    
    cpdef double distance(self, Coordinate p)
    cpdef bool is_nan(self)


cdef double radians(double degree)
cpdef tuple PLAP(Coordinate A, double L0, double a0, Coordinate B = *, bool inverse = *)
cpdef tuple PLLP(Coordinate A, double L0, double L1, Coordinate B, bool inverse = *)
cpdef tuple PLPP(Coordinate A, double L0, Coordinate B, Coordinate C, bool inverse = *)
cpdef tuple PXY(Coordinate A, double x, double y)

cdef bool legal_crank(Coordinate A, Coordinate B, Coordinate C, Coordinate D)
cdef str strbetween(str s, str front, str back)
cdef str strbefore(str s, str front)
