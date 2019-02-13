# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable classes of the validation in algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import (
    array as np_array,
    float32 as np_float32,
)
# Not a number and a large fitness. Infinity cannot be used for a chart.
from libc.math cimport NAN, HUGE_VAL
from numpy cimport ndarray
from verify cimport Verification
from expression cimport VPoint
from triangulation cimport vpoints_configure
from tinycadlib cimport (
    Coordinate,
    radians,
    plap,
    pllp,
    plpp,
    pxy,
    str_between,
    str_before,
)


@cython.final
cdef class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef ndarray upper, lower

    def __cinit__(self, dict mech_params):
        """mech_params = {
            'Expression': List[VPoint],
            'input': [(b0, d0), ...],
            'Placement': {pt: (x, y, r)},
            'Target': {pt: [(x0, y0), (x1, y1), ...]},
            'upper': ndarray[double, ndim=1],
            'lower': ndarray[double, ndim=1],
        }
        """
        # TODO: Get options.
        self.vpoints = list(mech_params.get('Expression', []))
        self.inputs = list(mech_params.get('input', []))
        self.exprs = vpoints_configure(self.vpoints, self.inputs)

        # Swap sorting.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]

    cdef ndarray[double, ndim=1] get_upper(self):
        return self.upper

    cdef ndarray[double, ndim=1] get_lower(self):
        return self.lower

    cdef int length(self):
        return len(self.upper)

    cdef double fitness(self, ndarray[double, ndim=1] v):
        """Chromosome format: (decided by upper and lower)
        
        v: [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        target: a list of target. [(1,5), (2,5), (3,5)]
        target_count: length of target.
        vars: mechanism variables count.
        """
        # TODO: fitness

    cpdef object result(self, ndarray[double, ndim=1] v):
        """Return the last answer."""
        # TODO: result
