# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Normalized planar four-bar linkage synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import array
from pyslvs.metaheuristics.utility cimport ObjFunc


cdef double trapezoidal_camp(double[:] a, double[:] b):
    """Error comparison by trapezoidal rule."""
    cdef double area = 0
    cdef int i
    for i in range(1, len(a)):
        area += abs(a[i - 1] + a[i] - b[i - 1] + b[i]) / 2
    return area


@cython.final
cdef class NPlanar(ObjFunc):
    """A normalized matching method.

    Defects free.
    """
    cdef double[:, :] target

    def __cinit__(self, dict mech):
        self.target = array(mech['target'])
        # TODO: Normalization

    cdef double fitness(self, double[:] v) nogil:
        pass

    cpdef object result(self, double[:] v):
        pass
