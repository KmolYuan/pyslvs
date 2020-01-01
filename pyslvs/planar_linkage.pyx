# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True, cdivision=True

"""The callable classes of the validation in the algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from collections import OrderedDict
from numpy import array as np_array, float64 as np_float
from libc.math cimport HUGE_VAL, NAN
from libcpp.list cimport list as clist
from .metaheuristics.utility cimport Objective
from .expression cimport get_vlinks, VJoint, VPoint, VLink
from .triangulation cimport (
    vpoints_configure,
    symbol_str,
    Expression,
    PLA,
    PLLP,
    PLPP,
    PXY,
)
from .bfgs cimport SolverSystem
from .tinycadlib cimport (
    radians,
    Coordinate,
    plap,
    pllp,
    plpp,
    pxy,
)


@cython.final
cdef class Planar(Objective):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef bint bfgs_mode
    cdef int target_count, v_base
    cdef clist[Expression] exprs
    cdef list vpoints, mapping_list
    cdef dict placement, target, mapping, mapping_r, data_dict
    cdef object inputs
    cdef double[:] upper, lower
    cdef SolverSystem bfgs_solver

    def __cinit__(self, dict mech):
        """mech = {
            'expression': List[VPoint],
            'input': OrderedDict([((b0, d0), [start, end]), ...]),
            'placement': {pt: (x, y, r)},
            'target': {pt: [(x0, y0), (x1, y1), ...]},
            'same': {pt: match_to_pt},
            # Bound has no position data.
            'upper': List[float],
            'lower': List[float],
        }
        """
        placement = mech.get('placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")
        target = mech.get('target', {})
        if len(target) == 0:
            raise ValueError("no target joint")
        check_set = set(map(len, target.values()))
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()
        # Change the target paths into memory view.
        self.target = {}
        same = mech.get('same', {})
        cdef int i, j
        cdef Coordinate[:] path
        for i in target:
            path = np_array([Coordinate(x, y) for x, y in target[i]], dtype=object)
            for j in range(i):
                if j in same:
                    i -= 1
            self.target[i] = path
        # Expressions
        self.vpoints = list(mech.get('expression', []))
        self.inputs = OrderedDict(mech.get('input', {}))
        status = {}
        self.exprs = vpoints_configure(self.vpoints, tuple(self.inputs.keys()), status).stack
        self.bfgs_mode = not all(status.values())
        # BFGS solver mode
        self.bfgs_solver = None
        # Data mapping
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_r = {v: k for k, v in self.mapping.items()}
        self.mapping_list = []
        # Bounds
        upper = list(mech.get('upper', []))
        lower = list(mech.get('lower', []))
        if len(upper) != len(lower):
            raise ValueError("upper and lower should be in the same size")
        # Position
        j = 0
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            upper[j:j] = (x + r, y + r)
            lower[j:j] = (x - r, y - r)
            self.mapping_list.append(i)
            j += 2
        self.v_base = len(upper)
        # Length of links
        cdef int a, b, c, d
        cdef VLink vlink
        for vlink in get_vlinks(self.vpoints):
            if len(vlink.points) < 2:
                continue
            if vlink.name == VLink.FRAME:
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
        # Input nodes and angle range
        for a, ((i, j), (start, end)) in enumerate(self.inputs.items()):
            upper.append(start)
            lower.append(end)
            for a in range(i):
                if a in same:
                    i -= 1
            for a in range(j):
                if a in same:
                    j -= 1
        upper[self.v_base:] *= self.target_count
        lower[self.v_base:] *= self.target_count
        self.upper = np_array(upper, dtype=np_float)
        self.lower = np_array(lower, dtype=np_float)
        # Swap upper and lower bound if reversed.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]
        # Result list
        self.data_dict = {}

    cpdef double[:] get_upper(self):
        return self.upper

    cpdef double[:] get_lower(self):
        return self.lower

    cpdef bint is_two_kernel(self):
        return self.bfgs_mode

    cdef inline double get_len(self, str expr1, str expr2):
        """Get the link length."""
        return self.mapping[frozenset({self.mapping_r[expr1], self.mapping_r[expr2]})]

    cdef inline bint solve(self, double[:] input_list):
        """Start solver function."""
        self.data_dict.clear()

        cdef int i
        cdef VPoint vpoint
        cdef Coordinate coord1, coord2
        for i, vpoint in enumerate(self.vpoints):
            if not vpoint.grounded():
                continue
            coord1 = Coordinate.__new__(Coordinate, vpoint.c[0][0], vpoint.c[0][1])
            if vpoint.type == VJoint.R:
                self.data_dict[self.mapping[i]] = coord1
                self.data_dict[i, -1] = coord1
            else:
                coord2 = Coordinate.__new__(Coordinate, vpoint.c[1][0], vpoint.c[1][1])
                self.data_dict[self.mapping[i]] = coord2
                self.data_dict[i, -1] = coord1
                self.data_dict[i, -2] = coord2

        # Solve
        i = 0
        cdef int t, params_count
        cdef Coordinate coord, coord3
        cdef Expression expr
        for expr in self.exprs:
            coord = Coordinate.__new__(Coordinate, NAN, NAN)
            if expr.func == PLA:
                target = symbol_str(expr.c2)
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord = plap(
                    coord1,
                    self.get_len(symbol_str(expr.c1), target),
                    radians(input_list[i])
                )
                i += 1
            elif expr.func == PLLP:
                target = symbol_str(expr.c3)
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord2 = self.data_dict[symbol_str(expr.c2)]
                coord = pllp(
                    coord1,
                    self.get_len(symbol_str(expr.c1), target),
                    self.get_len(symbol_str(expr.c2), target),
                    coord2,
                    expr.op
                )
            elif expr.func == PLPP:
                target = symbol_str(expr.c4)
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord2 = self.data_dict[symbol_str(expr.c2)]
                coord3 = self.data_dict[symbol_str(expr.c3)]
                coord = plpp(
                    coord1,
                    self.get_len(symbol_str(expr.c1), target),
                    coord2,
                    coord3,
                    expr.op
                )
            elif expr.func == PXY:
                target = symbol_str(expr.c2)
                vpoint = self.vpoints[self.mapping_r[target]]
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord = pxy(coord1, vpoint.c[0][0] - coord1.x, vpoint.c[0][1] - coord1.y)
            else:
                return False

            if coord.is_nan():
                return False

            t = self.mapping_r[target]
            vpoint = self.vpoints[t]
            self.data_dict[target] = coord
            if vpoint.type == VJoint.R:
                self.data_dict[t, -1] = coord
            else:
                self.data_dict[t, -1] = vpoint.c[0]
                self.data_dict[t, -2] = coord

        if not self.bfgs_mode:
            return True

        # Add coordinate of known points.
        p_data_dict = {}
        for i in range(len(self.vpoints)):
            # {1: 'A'} vs {'A': (10., 20.)}
            if self.mapping[i] in self.data_dict:
                p_data_dict[i] = self.data_dict[self.mapping[i]]

        # Add specified link lengths.
        for k, v in self.mapping.items():
            if type(k) is frozenset:
                p_data_dict[k] = v

        if self.bfgs_solver is None:
            self.bfgs_solver = SolverSystem(self.vpoints, {}, p_data_dict)
        else:
            self.bfgs_solver.set_data(p_data_dict)

        # Solve
        try:
            solved_bfgs = self.bfgs_solver.solve()
        except ValueError:
            return False

        # Format:
        # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
        # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
        for i in range(len(self.vpoints)):
            if self.mapping[i] in self.data_dict:
                continue
            vpoint = self.vpoints[i]
            # These points solved by Sketch Solve.
            if vpoint.type == VJoint.R:
                self.data_dict[i, -1] = Coordinate.__new__(
                    Coordinate,
                    solved_bfgs[i][0],
                    solved_bfgs[i][1]
                )
            else:
                self.data_dict[i, -1] = Coordinate.__new__(
                    Coordinate,
                    solved_bfgs[i][0][0],
                    solved_bfgs[i][0][1]
                )
                self.data_dict[i, -2] = Coordinate.__new__(
                    Coordinate,
                    solved_bfgs[i][1][0],
                    solved_bfgs[i][1][1]
                )

        return True

    cdef double fitness(self, double[:] v):
        """Chromosome format: (decided by upper and lower)

        v: [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        """
        cdef int target_index = 0
        cdef VPoint vpoint
        for m in self.mapping_list:
            if type(m) is int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[target_index], v[target_index + 1])
                target_index += 2
            else:
                self.mapping[m] = v[target_index]
                target_index += 1

        cdef double fitness = 0.

        cdef int index, node
        cdef Coordinate coord
        cdef Coordinate[:] path
        for target_index in range(self.target_count):
            index = self.v_base + target_index
            if not self.solve(v[index:index + self.target_count]):
                return HUGE_VAL
            for node, path in self.target.items():
                coord = self.data_dict[node, -1]
                fitness += coord.distance(path[target_index])

        return fitness

    cpdef object result(self, double[:] v):
        """Return the last answer."""
        cdef int target_index = 0
        cdef VPoint vpoint
        for m in self.mapping_list:
            if type(m) is int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[target_index], v[target_index + 1])
                target_index += 2
            else:
                self.mapping[m] = v[target_index]
                target_index += 1

        self.solve(v[self.v_base:self.v_base + self.target_count])
        expressions = []
        cdef int i
        cdef double x1, y1, x2, y2
        cdef Coordinate coord1, coord2
        for i in range(len(self.vpoints)):
            vpoint = self.vpoints[i]
            coord1 = self.data_dict[i, -1]
            vpoint.locate(coord1.x, coord1.y)
            if vpoint.type != VJoint.R:
                coord2 = self.data_dict[i, -2]
                vpoint.move((coord1.x, coord1.y), (coord2.x, coord2.y))
            expressions.append(vpoint.expr())

        return "M[" + ", ".join(expressions) + "]"
