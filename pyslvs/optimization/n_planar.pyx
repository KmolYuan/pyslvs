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
from libc.math cimport M_PI
from numpy import array
from numpy.fft import fft
from pyslvs.expression cimport VPoint
from pyslvs.tinycadlib cimport c_uniform_path, uniform_expr
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

    Defects free. Normalized parameters are $[L_0, L_2, L_3, L_4, \\alpha]$.

    ![pxy](img/uniform_four_bar.png)
    """
    cdef int len
    cdef double[:, :] target

    def __cinit__(self, dict mech):
        self.target = array(mech['target'])
        self.len = len(self.target)
        # TODO: Normalization
        self.lb = array([1e-5] * 4 + [0.])
        self.ub = array([5.] * 4 + [2. * M_PI])
        raise NotImplementedError

    cdef double fitness(self, double[:] v) nogil:
        """Generate linkage with 5 parameters."""
        cdef double[:, :] p = c_uniform_path(v[None, :], self.len)[0]
        # TODO: Normalization
        raise NotImplementedError

    cpdef object result(self, double[:] v):
        return "M[" + ", ".join([(<VPoint> vp).expr()
                                 for vp in uniform_expr(v)]) + "]"
