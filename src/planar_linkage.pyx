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
    sort as np_sort,
)
# Not a number and a large fitness. Infinity cannot be used for a chart.
from libc.math cimport hypot, HUGE_VAL
from numpy cimport ndarray
from verify cimport Verification
from expression cimport (
    get_vlinks,
    VPoint,
    VLink,
)
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
    expr_solving,
)


@cython.final
cdef class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef int target_count, base_index
    cdef list vpoints, inputs, exprs, mapping_list
    cdef dict placement, target, mapping
    cdef ndarray upper, lower

    def __cinit__(self, mech_params: dict):
        """mech_params = {
            'Expression': List[VPoint],
            'input': [(b0, d0), ...],
            'Placement': {pt: (x, y, r)},
            'Target': {pt: [(x0, y0), (x1, y1), ...]},
            # Bound has no position data.
            'upper': List[float],
            'lower': List[float],
        }
        """
        cdef dict placement = mech_params.get('Placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")

        self.target = mech_params.get('Target', {})
        if len(self.target) == 0:
            raise ValueError("no target joint")

        cdef set check_set = set(map(len, self.target.values()))
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()

        # Options
        self.vpoints = list(mech_params.get('Expression', []))
        self.inputs = list(mech_params.get('input', []))
        self.exprs = vpoints_configure(self.vpoints, self.inputs)

        # Bound
        cdef list upper = list(mech_params.get('upper', []))
        cdef list lower = list(mech_params.get('lower', []))
        if len(upper) != len(lower):
            raise ValueError("upper and lower should be in the same size")

        cdef list upper_input = []
        cdef list lower_input = []

        # Data mapping
        cdef int i
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_list = []

        # Position
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            upper_input.extend((x + r, y + r))
            lower_input.extend((x - r, y - r))
            self.mapping_list.append(i)

        # Length of links
        cdef int a, b, c, d
        cdef VLink vlink
        for vlink in get_vlinks(self.vpoints):
            if len(vlink.points) < 2:
                continue
            a = vlink.points[0]
            b = vlink.points[1]
            self.mapping[a, b] = None
            self.mapping_list.append((a, b))
            for c in vlink.points[2:]:
                for d in (a, b):
                    self.mapping[c, d] = None
                    self.mapping_list.append((c, d))

        upper_input.extend(upper[:-len(self.inputs)])
        lower_input.extend(lower[:-len(self.inputs)])

        for i in range(1, len(self.inputs) + 1):
            upper_input.extend([upper[-i]] * self.target_count)
            lower_input.extend([lower[-i]] * self.target_count)
        self.upper = np_array(upper_input, dtype=np_float32)
        self.lower = np_array(lower_input, dtype=np_float32)
        self.base_index = len(self.upper) - len(self.inputs) * self.target_count

        # Swap sorting.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]

    cdef ndarray[double, ndim=1] get_upper(self):
        return self.upper

    cdef ndarray[double, ndim=1] get_lower(self):
        return self.lower

    cdef double fitness(self, ndarray[double, ndim=1] v):
        """Chromosome format: (decided by upper and lower)
        
        v: [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        target: a list of target. [(1,5), (2,5), (3,5)]
        target_count: length of target.
        vars: mechanism variables count.
        """
        cdef int index = 0
        cdef VPoint vpoint
        cdef object m
        for m in self.mapping_list:
            if type(m) == int:
                vpoint = self.vpoints[m]
                vpoint.move((v[index], v[index + 1]))
                index += 2
            else:
                self.mapping[m] = v[index]
                index += 1

        cdef double fitness = 0.
        cdef ndarray[double, ndim=1] input_list = np_sort(v[self.base_index:])

        cdef int node, target_index
        cdef double x, y, tx, ty
        cdef list path, result
        cdef object coord
        for target_index in range(len(self.target)):
            try:
                result = expr_solving(
                    self.exprs,
                    self.mapping,
                    self.vpoints,
                    input_list[target_index::self.target_count]
                )
            except ValueError:
                return HUGE_VAL

            for node, path in self.target.items():
                for tx, ty in path:
                    coord = result[node]
                    if type(coord[0]) == tuple:
                        x, y = coord[1]
                    else:
                        x, y = coord
                    fitness += hypot(x - tx, y - ty)

        return fitness

    cpdef object result(self, ndarray[double, ndim=1] v):
        """Return the last answer."""
        # TODO: result
        return {}
