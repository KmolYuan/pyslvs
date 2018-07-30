# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable class of the validation in algorithm."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

import numpy as np
cimport numpy as np


cdef enum limit:
    maxGen, minFit, maxTime


cdef class Chromosome:
    cdef public int n
    cdef public double f
    cdef public np.ndarray v
    
    cdef double distance(self, Chromosome)
    cpdef void assign(self, Chromosome)


cdef class Verification:
    cdef np.ndarray get_upper(self)
    cdef np.ndarray get_lower(self)
    cdef int get_nParm(self)
    cpdef dict get_coordinates(self, np.ndarray v)
