# -*- coding: utf-8 -*-
# cython: language_level=3

from numpy cimport int16_t
cdef int16_t[:, :] product(tuple pool, object stop_func)
