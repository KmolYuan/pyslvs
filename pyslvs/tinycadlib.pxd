# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from .topo_config cimport EStack

cdef extern from "swappable_pair.hpp" nogil:
    cppclass SwappablePair:
        int first, second
        bint operator==(const SwappablePair & rhs) const
        bint operator!=(const SwappablePair & rhs) const

cdef double radians(double degree) nogil
cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints,
                        object angles=*)
