# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable class of the validation in algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

from numpy cimport ndarray


cdef class Verification:
    cdef ndarray[double, ndim=1] get_upper(self)
    cdef ndarray[double, ndim=1] get_lower(self)
    cdef double fitness(self, ndarray[double, ndim=1] v)
    cpdef object result(self, ndarray[double, ndim=1] v)
