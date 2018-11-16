# -*- coding: utf-8 -*-
# cython: language_level=3

"""Planarity check function for Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

from cpython cimport bool
from graph cimport Graph


cpdef bool is_planar(Graph g)
