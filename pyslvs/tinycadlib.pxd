# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from .topo_config cimport EStack

ctypedef (double, double) CCoord

cdef double radians(double degree) nogil
cdef CCoord cpxy(CCoord c1, double x, double y) nogil
cdef CCoord cppp(CCoord c1, CCoord c2, CCoord c3) nogil
cdef CCoord cplap(CCoord c1, double d0, double a0, CCoord c2,
                  bint inverse) nogil
cdef CCoord cpllp(CCoord c1, double d0, double d1, CCoord c2,
                  bint inverse) nogil
cdef CCoord cplpp(CCoord c1, double d0, CCoord c2, CCoord c3,
                  bint inverse) nogil
cdef CCoord cpalp(CCoord c1, double a0, double d0, CCoord c2,
                  bint inverse) nogil
cpdef void expr_parser(EStack exprs, dict data_dict)
cpdef tuple data_collecting(EStack exprs, dict mapping, object vpoints_)
cpdef list expr_solving(EStack exprs, dict mapping, object vpoints,
                        object angles=*)
