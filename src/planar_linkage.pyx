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
    zeros as np_zeros,
    array as np_array,
    float as np_float,
    sort as np_sort,
)
# Not a number and a large fitness. Infinity cannot be used for a chart.
from libc.math cimport (
    hypot,
    HUGE_VAL,
    isnan,
)
from numpy cimport ndarray
from verify cimport Verification
from expression cimport (
    get_vlinks,
    VJoint,
    VPoint,
    VLink,
)
from triangulation cimport vpoints_configure
from bfgs cimport vpoint_solving
from tinycadlib cimport (
    radians,
    expr_parser,
    data_collecting,
)


@cython.final
cdef class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef int target_count, base_index
    cdef list vpoints, inputs, exprs, mapping_list
    cdef dict placement, target, mapping
    cdef ndarray upper, lower, result_list

    def __cinit__(self, mech_params: dict):
        """mech_params = {
            'Expression': List[VPoint],
            'input': [(b0, d0), ...],
            'Placement': {pt: (x, y, r)},
            'Target': {pt: [(x0, y0), (x1, y1), ...]},
            'same': {pt: match_to_pt},
            # Bound has no position data.
            'upper': List[float],
            'lower': List[float],
        }
        """
        cdef dict placement = mech_params.get('Placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")

        cdef dict target = mech_params.get('Target', {})
        if len(target) == 0:
            raise ValueError("no target joint")

        cdef set check_set = set(map(len, target.values()))
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()

        # Change the target paths into memory view.
        self.target = {}
        cdef dict same = mech_params.get('same', {})

        cdef int i, j
        cdef double[:, :] path
        for i in target:
            path = np_array(target[i], dtype=np_float)
            for j in range(i):
                if j in same:
                    i -= 1
            self.target[i] = path

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
        self.upper = np_array(upper_input, dtype=np_float)
        self.lower = np_array(lower_input, dtype=np_float)
        self.base_index = len(self.upper) - len(self.inputs) * self.target_count

        # Swap sorting.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]
        
        # Result list
        self.result_list = np_zeros((len(self.vpoints), 2, 2), dtype=np_float)

    cdef ndarray[double, ndim=1] get_upper(self):
        return self.upper

    cdef ndarray[double, ndim=1] get_lower(self):
        return self.lower

    cdef inline bint solve(self, double[:] input_list, bint no_slider):
        """Start solver function."""
        # TODO: Need to be optimized.
        cdef dict data_dict
        cdef int dof
        data_dict, dof = data_collecting(self.exprs, self.mapping, self.vpoints)

        # Angles.
        cdef double a
        cdef int i
        for i, a in enumerate(input_list):
            data_dict[f'a{i}'] = radians(a)

        # Solve
        expr_parser(self.exprs, data_dict)

        cdef dict p_data_dict = {}
        cdef bint has_not_solved = False

        # Add coordinate of known points.
        for i in range(len(self.vpoints)):
            # {1: 'A'} vs {'A': (10., 20.)}
            if self.mapping[i] in data_dict:
                p_data_dict[i] = data_dict[self.mapping[i]]
            else:
                has_not_solved = True

        # Calling Sketch Solve kernel and try to get the result.
        cdef list solved_bfgs = []
        if has_not_solved:

            # Add specified link lengths.
            for k, v in data_dict.items():
                if type(k) == tuple:
                    p_data_dict[k] = v

            # Solve
            try:
                solved_bfgs = vpoint_solving(self.vpoints, {}, p_data_dict)
            except ValueError:
                return False

        # Format:
        # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
        # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
        cdef VPoint vpoint
        for i in range(len(self.vpoints)):
            vpoint = self.vpoints[i]
            if self.mapping[i] in data_dict:
                # These points has been solved.
                if isnan(data_dict[self.mapping[i]][0]):
                    return False
                if no_slider or vpoint.type == VJoint.R:
                    self.result_list[i, 0] = data_dict[self.mapping[i]]
                else:
                    self.result_list[i, 0] = vpoint.c[0]
                    self.result_list[i, 1] = data_dict[self.mapping[i]]
            elif solved_bfgs:
                # These points solved by Sketch Solve.
                if no_slider or vpoint.type == VJoint.R:
                    self.result_list[i, 0] = solved_bfgs[i]
                else:
                    self.result_list[i, 0] = solved_bfgs[i][0]
                    self.result_list[i, 1] = solved_bfgs[i][1]
            else:
                # No answer.
                if no_slider or vpoint.type == VJoint.R:
                    self.result_list[i, 0] = vpoint.c[0]
                else:
                    self.result_list[i, 0] = vpoint.c[0]
                    self.result_list[i, 1] = vpoint.c[1]

        return True

    cdef double fitness(self, ndarray[double, ndim=1] v):
        """Chromosome format: (decided by upper and lower)

        v: [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        """
        cdef int target_index = 0
        cdef VPoint vpoint
        cdef object m
        for m in self.mapping_list:
            if type(m) == int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[target_index], v[target_index + 1])
                target_index += 2
            else:
                self.mapping[m] = v[target_index]
                target_index += 1

        cdef double[:] input_list = np_sort(v[self.base_index:])
        cdef double fitness = 0.

        cdef int node
        cdef double x, y, tx, ty
        cdef double[:, :] path
        for target_index in range(self.target_count):
            if not self.solve(input_list[target_index::self.target_count], True):
                return HUGE_VAL

            for node, path in self.target.items():
                tx, ty = path[target_index]
                x, y = self.result_list[node, 0]
                fitness += hypot(x - tx, y - ty)

        return fitness

    cpdef object result(self, ndarray[double, ndim=1] v):
        """Return the last answer."""
        cdef int target_index = 0
        cdef VPoint vpoint
        cdef object m
        for m in self.mapping_list:
            if type(m) == int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[target_index], v[target_index + 1])
                target_index += 2
            else:
                self.mapping[m] = v[target_index]
                target_index += 1

        cdef double[:] input_list = np_sort(v[self.base_index:])
        self.solve(input_list[::self.target_count], False)

        cdef list expressions = []

        cdef int i
        cdef double x1, y1, x2, y2
        for i in range(len(self.vpoints)):
            vpoint = self.vpoints[i]
            x1, y1 = self.result_list[i, 0]
            vpoint.locate(x1, y1)
            if vpoint.type != VJoint.R:
                x1, y1 = self.result_list[i, 0]
                x2, y2 = self.result_list[i, 0]
                vpoint.move((x1, y1), (x2, y2))
            expressions.append(vpoint.expr())

        return "M[" + ", ".join(expressions) + "]"
