# -*- coding: utf-8 -*-
# cython: language_level=3

"""PMKS simbolics."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from cpython cimport bool
from numpy cimport ndarray

cdef class VPoint:
    
    # VPoint(links, type_int, angle, color_str, x, y, color_func=None)
    
    # Members
    cdef readonly tuple links
    cdef readonly ndarray c
    cdef readonly int type
    cdef readonly object color
    cdef readonly str colorSTR
    cdef readonly str typeSTR
    cdef readonly double x, y, angle
    cdef double __offset
    cdef bool __has_offset
    
    # Set values
    cpdef void move(self, tuple, tuple c2 = *) except *
    cpdef void rotate(self, double)
    cpdef void set_offset(self, double)
    
    # Get or calculate values.
    cpdef double distance(self, VPoint)
    cpdef bool has_offset(self)
    cpdef double offset(self)
    cpdef double true_offset(self)
    cpdef double slope_angle(self, VPoint, int num1 = *, int num2 = *)
    
    # Link operators.
    cpdef bool grounded(self)
    cpdef bool pin_grounded(self)
    cpdef bool same_link(self, VPoint p)
    cpdef bool no_link(self)
    cpdef bool is_slot_link(self, str)
