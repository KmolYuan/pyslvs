# -*- coding: utf-8 -*-
# cython: language_level=3

"""Number synthesis.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from numpy cimport int16_t

cdef int16_t[:, :] product(tuple pool, object stop_func)
