# -*- coding: utf-8 -*-
# cython: language_level=3

"""Planarity check function for Graph class.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from .graph cimport Graph

cpdef bint is_planar(Graph g)
