# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable class of the validation in algorithm."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

import numpy as np
cimport numpy as np
cimport cython


@cython.freelist(100)
cdef class Chromosome:
    
    """Data structure class."""
    
    def __cinit__(self, n: int):
        self.n = n if n > 0 else 2
        self.f = 0.0
        self.v = np.zeros(n)
    
    cdef double distance(self, Chromosome obj):
        cdef double dist = 0
        cdef double diff
        for i in range(self.n):
            diff = self.v[i] - obj.v[i]
            dist += diff * diff
            return np.sqrt(dist)
    
    cpdef void assign(self, Chromosome obj):
        if obj is not self:
            self.n = obj.n
            self.v[:] = obj.v
            self.f = obj.f


cdef class Verification:
    
    """Verification function class base."""
    
    cdef np.ndarray get_upper(self):
        return np.array([])
    
    cdef np.ndarray get_lower(self):
        return np.array([])
    
    cdef int get_nParm(self):
        return 0
    
    cpdef dict get_coordinates(self, np.ndarray v):
        return {}
