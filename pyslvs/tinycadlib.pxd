# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.vector cimport vector
from libcpp.map cimport map
from .topo_config cimport EStack
from .expression cimport CCoord

cdef extern from "swappable_pair.hpp" nogil:
    cppclass SwappablePair:
        int first, second
        bint operator ==(const SwappablePair & rhs) const
        bint operator !=(const SwappablePair & rhs) const

cdef double radians(double degree) nogil
cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints,
                        object angles=*)

cdef class ExprSolver:
    cdef EStack exprs
    cdef vector[CCoord] joint_pos
    cdef map[SwappablePair, double] link_length
    # TODO: STL container

    cdef double[:, :] solve(self) nogil
