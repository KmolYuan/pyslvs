# -*- coding: utf-8 -*-
# cython: language_level=3

"""Wrapper header of BFGS algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.list cimport list as clist
from libcpp.vector cimport vector
from libcpp.map cimport map
from .sketch_solve cimport Point, Line, Constraint
from .expression cimport VLink

cdef class SolverSystem:
    # SolverSystem(object vpoints, dict inputs = None, dict data_dict = None)

    cdef list vpoints
    cdef dict vlinks
    cdef dict inputs
    cdef dict data_dict

    cdef clist[double] params
    cdef clist[double] constants
    cdef vector[Point] points
    cdef map[int, int] sliders
    cdef vector[Point] slider_bases
    cdef vector[Point] slider_slots
    cdef clist[Line] slider_lines

    cdef clist[Constraint] cons_list
    cdef clist[Line] handles
    cdef clist[double] inputs_angle  # For editing angles
    cdef clist[double] data_values  # For editing custom values

    cpdef bint same_points(self, object vpoints_)
    cpdef frozenset show_inputs(self)
    cpdef frozenset show_data(self)
    cdef Point *point_ptr(self, int i, VLink vlink)
    cdef void build_expression(self)
    cdef bint check_known(self, int i)
    cpdef void set_inputs(self, dict inputs)
    cpdef void set_data(self, dict data_dict)
    cpdef list solve(self)
