# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True, cdivision=True

"""Differential Evolution

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import array, float64 as np_float
from pyslvs.metaheuristics.verify cimport Verification

@cython.final
cdef class TestObj(Verification):

    """Test objective function."""

    def __cinit__(self):
        pass

    cdef double[:] get_upper(self):
        return array([100, 100], dtype=np_float)

    cdef double[:] get_lower(self):
        return array([0, 0], dtype=np_float)

    cdef double fitness(self, double[:] v):
        cdef double x1 = v[0]
        cdef double x2 = v[1]
        return x1 * x1 + 8 * x2

    cpdef object result(self, double[:] v):
        cdef double x1 = v[0]
        cdef double x2 = v[1]
        return tuple(v), x1 * x1 + 8 * x2
