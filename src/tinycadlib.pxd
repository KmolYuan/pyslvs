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
    
    cpdef double distance(self, Coordinate)
    cpdef bool is_nan(self)


cpdef tuple PLAP(Coordinate, double, double, Coordinate B = *, bool inverse = *)
cpdef tuple PLLP(Coordinate, double, double, Coordinate, bool inverse = *)
cpdef tuple PLPP(Coordinate, double, Coordinate, Coordinate, bool inverse = *)
cpdef tuple PXY(Coordinate, double, double)

cdef bool legal_crank(Coordinate, Coordinate, Coordinate, Coordinate)
cdef str strbetween(str, str, str)
cdef str strbefore(str, str)
