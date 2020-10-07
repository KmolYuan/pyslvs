# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from .topo_config cimport EStack

cdef double radians(double degree) nogil
cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints,
                        object angles=*)

cdef class ExprSolver:
    cdef EStack exprs
    # TODO: STL container

    cdef double[:, :] solve(self) nogil
