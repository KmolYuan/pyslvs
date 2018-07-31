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
    """def __cinit__(self,
        links: str,
        type_int: int,
        angle: double,
        color_str: str,
        x: double,
        y: double,
        color_func: object = None
    )
    """
    cdef readonly tuple links
    cdef readonly ndarray c
    cdef readonly int type
    cdef readonly object color
    cdef readonly str colorSTR
    cdef readonly str typeSTR
    cdef readonly double x, y, angle
    
    cpdef void move(self, tuple, tuple c2 = *)
    cpdef void rotate(self, double)
    cpdef double distance(self, VPoint)
    cpdef double slope_angle(self, VPoint, int num1 = *, int num2 = *)
    
    cpdef bool grounded(self)
    cpdef bool no_link(self)
    cpdef bool is_slot_link(self, str)
