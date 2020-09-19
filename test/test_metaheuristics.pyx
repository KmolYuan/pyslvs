# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False

"""Differential Evolution

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import array, float64 as np_float
from unittest import TestCase
from pyslvs.metaheuristics import ALGORITHM, PARAMS
from pyslvs.metaheuristics.utility cimport ObjFunc


@cython.final
cdef class TestObj(ObjFunc):
    """Test objective function.

    f(x) = x1^2 + 8*x2
    """

    def __cinit__(self):
        self.ub = array([100, 100], dtype=np_float)
        self.lb = array([0, 0], dtype=np_float)

    cdef double target(self, double[:] v) nogil:
        return v[0] * v[0] + 8 * v[1]

    cdef double fitness(self, double[:] v) nogil:
        return self.target(v)

    cpdef object result(self, double[:] v):
        return tuple(v), self.target(v)


class AlgorithmTest(TestCase):

    def test_obj_func(self):
        """Test with an objective function."""
        settings = {'min_fit': 1e-20, 'report': 10}
        obj = TestObj()
        for t, setting in PARAMS.items():
            settings.update(setting)
            x, fval = ALGORITHM[t](obj, settings).run()
            self.assertAlmostEqual(0., x[0], 6)
            self.assertAlmostEqual(0., x[1], 6)
            self.assertAlmostEqual(0., fval, 6)
