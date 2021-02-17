# -*- coding: utf-8 -*-
# cython: language_level=3

"""Fast planar linkage synthesis. (defect allowable)

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cdef void roll(double[:, :] path, int ind) nogil
