# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True

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
    NAN,
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
    Coordinate,
    plap,
    pllp,
    plpp,
    pxy,
)


@cython.final
cdef class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef bint bfgs_mode
    cdef int target_count, base_index
    cdef tuple exprs
    cdef list vpoints, inputs, mapping_list
    cdef dict placement, target, mapping, mapping_r
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

        cdef list inputs = list(mech_params.get('input', []))

        # Change the target paths into memory view.
        self.target = {}
        self.inputs = []
        cdef dict same = mech_params.get('same', {})

        cdef int i, j
        cdef double[:, :] path
        for i in target:
            path = np_array(target[i], dtype=np_float)
            for j in range(i):
                if j in same:
                    i -= 1
            self.target[i] = path

        cdef int a
        for i, j in inputs:
            for a in range(i):
                if a in same:
                    i -= 1
            for a in range(j):
                if a in same:
                    j -= 1
            self.inputs.append((i, j))

        # Options
        self.vpoints = list(mech_params.get('Expression', []))
        cdef dict status = {}
        self.exprs = tuple(vpoints_configure(self.vpoints, self.inputs, status))
        self.bfgs_mode = not all(status.values())

        # Bound
        cdef list upper = list(mech_params.get('upper', []))
        cdef list lower = list(mech_params.get('lower', []))
        if len(upper) != len(lower):
            raise ValueError("upper and lower should be in the same size")

        cdef list upper_input = []
        cdef list lower_input = []

        # Data mapping
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_r = {v: k for k, v in self.mapping.items()}
        self.mapping_list = []

        # Position
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            upper_input.extend((x + r, y + r))
            lower_input.extend((x - r, y - r))
            self.mapping_list.append(i)

        # Length of links
        cdef int b, c, d
        cdef frozenset pair
        cdef VLink vlink
        for vlink in get_vlinks(self.vpoints):
            if len(vlink.points) < 2:
                continue
            a = vlink.points[0]
            b = vlink.points[1]
            pair = frozenset({a, b})
            self.mapping[pair] = None
            self.mapping_list.append(pair)
            for c in vlink.points[2:]:
                for d in (a, b):
                    pair = frozenset({c, d})
                    self.mapping[pair] = None
                    self.mapping_list.append(pair)

        upper_input.extend(upper[:-len(self.inputs)])
        lower_input.extend(lower[:-len(self.inputs)])

        for i in range(1, len(self.inputs) + 1):
            upper_input.extend([upper[-i]] * self.target_count)
            lower_input.extend([lower[-i]] * self.target_count)
        self.upper = np_array(upper_input, dtype=np_float)
        self.lower = np_array(lower_input, dtype=np_float)
        self.base_index = len(self.upper) - len(self.inputs) * self.target_count

        # Swap upper and lower bound if reversed.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]

        # Result list
        self.result_list = np_zeros((len(self.vpoints), 2, 2), dtype=np_float)

    cdef ndarray[double, ndim=1] get_upper(self):
        return self.upper

    cdef ndarray[double, ndim=1] get_lower(self):
        return self.lower

    cpdef bint is_two_kernel(self):
        return self.bfgs_mode

    cdef inline double get_len(self, str expr1, str expr2):
        """Get the link length."""
        return self.mapping[frozenset({self.mapping_r[expr1], self.mapping_r[expr2]})]

    cdef inline bint solve(self, double[:] input_list):
        """Start solver function."""
        # TODO: Need to be optimized.
        cdef dict data_dict = {}

        cdef int i
        cdef VPoint vpoint
        for i, vpoint in enumerate(self.vpoints):
            if not vpoint.grounded():
                continue
            if vpoint.type == VJoint.R:
                self.result_list[i, 0] = vpoint.c[0]
                data_dict[self.mapping[i]] = vpoint.c[0]
            else:
                self.result_list[i, 0] = vpoint.c[0]
                self.result_list[i, 1] = vpoint.c[1]
                data_dict[self.mapping[i]] = vpoint.c[1]

        # Solve
        cdef int t, params_count
        cdef double x, y, x1, y1, x2, y2, x3, y3, a, d1, d2
        cdef str func, target
        cdef tuple expr
        i = 0
        for expr in self.exprs:
            # If the mechanism has no any solution.
            if not expr:
                break

            func = expr[0]
            target = expr[-1]
            params_count = len(expr) - 2
            t = self.mapping_r[target]
            vpoint = self.vpoints[t]
            x = NAN
            y = NAN
            if func == 'PLAP':
                x1, y1 = data_dict[expr[1]]
                d1 = self.get_len(expr[1], target)
                a = radians(input_list[i])
                if params_count == 3:
                    x, y = plap(Coordinate(x1, y1), d1, a)
                else:
                    x2, y2 = data_dict[expr[4]]
                    if params_count == 4:
                        x, y = plap(Coordinate(x1, y1), d1, a, Coordinate(x2, y2))
                    elif params_count == 5:
                        x, y = plap(
                            Coordinate(x1, y1),
                            d1,
                            a,
                            Coordinate(x2, y2),
                            expr[5] == 'T'
                        )
                i += 1
            elif func == 'PLLP':
                x1, y1 = data_dict[expr[1]]
                d1 = self.get_len(expr[1], target)
                d2 = self.get_len(expr[4], target)
                x2, y2 = data_dict[expr[4]]
                if params_count == 4:
                    x, y = pllp(Coordinate(x1, y1), d1, d2, Coordinate(x2, y2))
                elif params_count == 5:
                    x, y = pllp(
                        Coordinate(x1, y1),
                        d1,
                        d2,
                        Coordinate(x2, y2),
                        expr[5] == 'T'
                    )
            elif func == 'PLPP':
                x1, y1 = data_dict[expr[1]]
                d1 = self.get_len(expr[1], target)
                x2, y2 = data_dict[expr[3]]
                x3, y3 = data_dict[expr[4]]
                if params_count == 4:
                    x, y = plpp(
                        Coordinate(x1, y1),
                        d1,
                        Coordinate(x2, y2),
                        Coordinate(x3, y3)
                    )
                elif params_count == 5:
                    x, y = plpp(
                        Coordinate(x1, y1),
                        d1,
                        Coordinate(x2, y2),
                        Coordinate(x3, y3),
                        expr[5] == 'T'
                    )
            elif func == 'PXY':
                x1, y1 = data_dict[expr[1]]
                x, y = pxy(Coordinate(x1, y1), vpoint.c[0][0] - x1, vpoint.c[0][1] - y1)

            if isnan(x):
                return False

            data_dict[target] = (x, y)
            if vpoint.type == VJoint.R:
                self.result_list[t, 0] = (x, y)
            else:
                self.result_list[t, 0] = vpoint.c[0]
                self.result_list[t, 1] = (x, y)

        if not self.bfgs_mode:
            return True

        # Calling Sketch Solve kernel and try to get the result.
        cdef dict p_data_dict
        # Add coordinate of known points.
        p_data_dict = {}
        for i in range(len(self.vpoints)):
            # {1: 'A'} vs {'A': (10., 20.)}
            if self.mapping[i] in data_dict:
                p_data_dict[i] = data_dict[self.mapping[i]]

        # Add specified link lengths.
        for k, v in self.mapping.items():
            if type(k) == frozenset:
                p_data_dict[k] = v

        # Solve
        cdef list solved_bfgs
        try:
            solved_bfgs = vpoint_solving(self.vpoints, {}, p_data_dict)
        except ValueError:
            return False

        # Format:
        # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
        # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
        for i in range(len(self.vpoints)):
            if self.mapping[i] in data_dict:
                continue

            vpoint = self.vpoints[i]
            # These points solved by Sketch Solve.
            if vpoint.type == VJoint.R:
                self.result_list[i, 0] = solved_bfgs[i]
            else:
                self.result_list[i, 0] = solved_bfgs[i][0]
                self.result_list[i, 1] = solved_bfgs[i][1]

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
            if not self.solve(input_list[target_index::self.target_count]):
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
        self.solve(input_list[::self.target_count])

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
