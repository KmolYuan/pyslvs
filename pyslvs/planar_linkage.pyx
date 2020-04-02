# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable classes of the validation in the algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from collections import OrderedDict
from numpy cimport ndarray
from numpy import zeros, array, float64 as np_float
from pywt import dwt
from libc.math cimport HUGE_VAL, M_PI, cos, sin, atan2, INFINITY as INF
from .metaheuristics.utility cimport Objective
from .expression cimport VJoint, VPoint
from .triangulation cimport (t_config, symbol_str, I_LABEL, A_LABEL, Expr,
    PLA, PLAP, PLLP, PLPP, PXY, EStack)
from .bfgs cimport SolverSystem
from .tinycadlib cimport radians, Coordinate, plap, pllp, plpp, pxy

DEF WAVELET = "db3"


def norm_path(path, scale=1):
    """Python wrapper of normalization function."""
    cdef Coordinate[:] path_m = array([
        Coordinate.__new__(Coordinate, x, y) for x, y in path], dtype=object)
    _normalization(path_m, scale)
    return [(c.x, c.y) for c in path_m]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _normalization(Coordinate[:] path, double scale):
    """Normalization implementation."""
    cdef Coordinate centre = Coordinate.__new__(Coordinate, 0, 0)
    cdef Coordinate c
    for c in path:
        centre.x += c.x
        centre.y += c.y
    cdef size_t end = len(path)
    centre.x /= end
    centre.y /= end
    cdef double[:] angle = zeros(end + 1, dtype=np_float)
    cdef double[:] length = zeros(end + 1, dtype=np_float)
    cdef size_t sp = 0
    cdef size_t i
    for i in range(len(path)):
        c = path[i]
        angle[i] = atan2(c.y - centre.y, c.x - centre.x)
        length[i] = centre.distance(c)
        if length[i] > length[end]:
            length[end] = length[i]
            angle[end] = angle[i]
            sp = end
    _aligned(path, sp)
    cdef double[:] bound = array([INF, -INF], dtype=np_float)
    cdef double a
    for i in range(len(path)):
        c = path[i]
        a = angle[i] - angle[end]
        c.x = length[i] * cos(a)
        c.y = length[i] * sin(a)
        if c.x < bound[0]:
            bound[0] = c.x
        if c.x > bound[1]:
            bound[1] = c.x
    scale /= (bound[1] - bound[0])
    for c in path:
        c.x *= scale
        c.y *= scale


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _aligned(Coordinate[:] path, size_t sp):
    """Split 1D path from sp, concatenate to end."""
    if sp == 0:
        return
    cdef double[:, :] tmp = zeros((sp, 2), dtype=np_float)
    cdef size_t i
    for i in range(sp):
        tmp[i, 0] = path[i].x
        tmp[i, 1] = path[i].y
    for i in range(sp, len(path)):
        path[i - sp].x = path[i].x
        path[i - sp].y = path[i].y
    for i in range(sp):
        path[len(path) - sp + i].x = tmp[i, 0]
        path[len(path) - sp + i].y = tmp[i, 1]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _wavelet(Coordinate[:] path) except *:
    """Return DWT."""
    cdef double[:, :] wave = zeros((len(path), 2), dtype=np_float)
    cdef size_t i
    cdef Coordinate c
    for i in range(len(path)):
        c = path[i]
        wave[i, 0] = c.x
        wave[i, 1] = c.y
    xa, xd = dwt(wave[:, 0], WAVELET)
    ya, yd = dwt(wave[:, 1], WAVELET)
    wave = zeros((4, max(len(xa), len(xd))), dtype=np_float)
    for i in range(len(xa)):
        wave[0, i] = xa[i]
        wave[1, i] = ya[i]
    for i in range(len(xd)):
        wave[2, i] = xd[i]
        wave[3, i] = yd[i]
    return wave


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _cmp_wavelet(double[:, :] wave1, double[:, :] wave2):
    """Compare two waves."""
    cdef double fitness = 0
    cdef size_t i
    for i in range(len(wave1)):
        fitness += (
            abs(wave1[0, i] - wave2[0, i])
            + abs(wave1[1, i] - wave2[1, i])
            + abs(wave1[2, i] - wave2[2, i]) * 1e3
            + abs(wave1[3, i] - wave2[3, i]) * 1e3
        )
    return fitness


@cython.final
@cython.boundscheck(False)
@cython.wraparound(False)
cdef class Planar(Objective):
    """This class is used to verified kinematics of the linkage mechanism."""
    cdef bint bfgs_mode, shape_only, wavelet_mode, ordered
    cdef int target_count, input_count, l_base
    cdef list vpoints, mapping_list
    cdef dict placement, target, mapping, mapping_r, data_dict
    cdef object inputs
    cdef double[:] upper, lower, polar_angles
    cdef EStack exprs
    cdef SolverSystem bfgs_solver

    def __cinit__(self, dict mech):
        # mech = {
        #     'expression': List[VPoint],
        #     'input': OrderedDict([((b0, d0), [start, end]), ...]),
        #     'placement': {pt: (x, y, r)},
        #     'target': {pt: [(x0, y0), (x1, y1), ...]},
        #     'same': {pt: match_to_pt},
        #     # Bounds of base link length
        #     'upper': float,
        #     'lower': float,
        #     'shape_only': bool,
        #     'wavelet_mode': bool,
        # }
        placement = mech.get('placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")
        target = mech.get('target', {})
        if len(target) == 0:
            raise ValueError("no target joint")
        check_set = {len(t) for t in target.values()}
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()
        # Change the target paths into memory view
        self.target = {}
        same = mech.get('same', {})
        self.shape_only = mech.get('shape_only', False)
        self.wavelet_mode = mech.get('wavelet_mode', False)
        cdef int i, j
        cdef double[:, :] wave
        cdef Coordinate[:] path
        for i in target:
            path = array([Coordinate(x, y) for x, y in target[i]], dtype=object)
            for j in range(i):
                if j in same:
                    i -= 1
            if self.shape_only or self.wavelet_mode:
                _normalization(path, 1)
                if self.wavelet_mode:
                    self.target[i] = _wavelet(path)
                    continue
            self.target[i] = path
        # Expressions
        self.vpoints = list(mech.get('expression', []))
        self.inputs = OrderedDict(mech.get('input', {}))
        status = {}
        self.exprs = t_config(self.vpoints, tuple(self.inputs.keys()), status)
        self.bfgs_mode = not all(status.values())
        # BFGS solver mode
        self.bfgs_solver = None
        # Data mapping
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_r = {v: k for k, v in self.mapping.items()}
        self.mapping_list = []
        # Bounds
        upper = []
        lower = []
        if len(upper) != len(lower):
            raise ValueError("upper and lower should be in the same size")
        # Position
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            upper.append(x + r)
            upper.append(y + r)
            lower.append(x - r)
            lower.append(y - r)
            self.mapping_list.append(i)
        cdef int p_base = len(upper)
        # Length of links
        link_upper = float(mech.get('upper', 100))
        link_lower = float(mech.get('lower', 0))
        i = 0
        cdef Expr expr
        for expr in self.exprs.stack:
            upper.append(link_upper)
            lower.append(link_lower)
            sym = frozenset({expr.c1.second, expr.target.second})
            self.mapping[sym] = None
            self.mapping_list.append(sym)
            if expr.func in {PLA, PLAP}:
                if expr.v2.first == A_LABEL:
                    upper.append(2 * M_PI)
                    lower.append(0)
                    sym = f"A{i}"
                    self.mapping[sym] = None
                    self.mapping_list.append(sym)
                    i += 1
            elif expr.func == PLLP:
                upper.append(link_upper)
                lower.append(link_lower)
                sym = frozenset({expr.c2.second, expr.target.second})
                self.mapping[sym] = None
                self.mapping_list.append(sym)
        self.l_base = len(upper)
        # Input nodes
        self.input_count = len(self.inputs)
        cdef int a
        for start, end in self.inputs.values():
            upper.append(radians(start))
            lower.append(radians(end))
        # Angle rage
        upper[self.l_base:] *= self.target_count
        lower[self.l_base:] *= self.target_count
        self.upper = array(upper, dtype=np_float)
        self.lower = array(lower, dtype=np_float)
        # Swap upper and lower bound if reversed
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]
        # Allocate memory
        self.polar_angles = zeros(self.l_base - p_base, dtype=np_float)
        # Result list
        self.data_dict = {}

    cpdef double[:] get_upper(self):
        """Return upper bound."""
        return self.upper

    cpdef double[:] get_lower(self):
        """Return lower bound."""
        return self.lower

    cpdef bint is_two_kernel(self):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        return self.bfgs_mode

    cdef inline double get_len(self, str expr1, str expr2):
        return self.mapping[frozenset({self.mapping_r[expr1], self.mapping_r[expr2]})]

    cdef bint solve(self, double[:] input_list):
        self.data_dict.clear()
        cdef int i
        cdef VPoint vpoint
        cdef Coordinate coord1, coord2, coord3
        for i, vpoint in enumerate(self.vpoints):
            coord1 = Coordinate.__new__(Coordinate, vpoint.c[0, 0], vpoint.c[0, 1])
            if vpoint.type == VJoint.R:
                self.data_dict[self.mapping[i]] = coord1
                self.data_dict[i, -1] = coord1
            else:
                coord2 = Coordinate.__new__(Coordinate, vpoint.c[1, 0], vpoint.c[1, 1])
                self.data_dict[self.mapping[i]] = coord2
                self.data_dict[i, -1] = coord1
                self.data_dict[i, -2] = coord2
        # Solve
        i = 0
        cdef int a = 0
        cdef int t, params_count
        cdef double length, angle
        cdef Expr expr
        for expr in self.exprs.stack:
            target = symbol_str(expr.target)
            if expr.func in {PLA, PLAP}:
                coord1 = self.data_dict[symbol_str(expr.c1)]
                if expr.func == PLAP:
                    coord2 = self.data_dict[symbol_str(expr.c2)]
                length = self.get_len(symbol_str(expr.c1), target)
                if expr.v2.first == I_LABEL:
                    angle = input_list[i]
                    i += 1
                else:
                    angle = self.polar_angles[a]
                    a += 1
                if expr.func == PLA:
                    coord3 = plap(coord1, length, angle)
                else:
                    coord3 = plap(coord1, length, angle, coord2, expr.op)
            elif expr.func == PLLP:
                coord3 = pllp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.get_len(symbol_str(expr.c1), target),
                    self.get_len(symbol_str(expr.c2), target),
                    self.data_dict[symbol_str(expr.c2)],
                    expr.op
                )
            elif expr.func == PLPP:
                coord3 = plpp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.get_len(symbol_str(expr.c1), target),
                    self.data_dict[symbol_str(expr.c2)],
                    self.data_dict[symbol_str(expr.c3)],
                    expr.op
                )
            elif expr.func == PXY:
                vpoint = self.vpoints[self.mapping_r[target]]
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord3 = pxy(
                    coord1,
                    vpoint.c[0, 0] - coord1.x,
                    vpoint.c[0, 1] - coord1.y
                )
            else:
                return False
            if coord3.is_nan():
                return False
            t = self.mapping_r[target]
            vpoint = self.vpoints[t]
            self.data_dict[target] = coord3
            if vpoint.type == VJoint.R:
                self.data_dict[t, -1] = coord3
            else:
                self.data_dict[t, -1] = (vpoint.c[0, 0], vpoint.c[0, 1])
                self.data_dict[t, -2] = coord3
        if not self.bfgs_mode:
            return True
        # Add coordinate of known points
        p_data_dict = {}
        for i in range(len(self.vpoints)):
            # {1: 'A'} vs {'A': (10., 20.)}
            if self.mapping[i] in self.data_dict:
                p_data_dict[i] = self.data_dict[self.mapping[i]]
        # Add specified link lengths
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
            # These points solved by Sketch Solve
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
        """Return fitness from chromosome.

        + Coordinates of fixed pivots. [0:self.l_base]
            [(xn, yn), ...]
        + Length and the angles of the links. [self.l_base]
        + Angle respect to the target points.
        """
        cdef int index = 0
        cdef int a_index = 0
        for m in self.mapping_list:
            if type(m) is int:
                (<VPoint>self.vpoints[m]).locate(v[index], v[index + 1])
                index += 2
            elif type(m) is str and m.startswith('A'):
                self.polar_angles[a_index] = v[index]
                a_index += 1
                index += 1
            else:
                self.mapping[m] = v[index]
                index += 1
        cdef ndarray[double, ndim=2] angles = zeros(
            (self.input_count, self.target_count), dtype=np_float)
        for index in range(self.input_count):
            a_index = index + self.l_base
            angles[index, :] = array(v[a_index:a_index + self.target_count])
        if self.wavelet_mode:
            angles.sort(axis=1)
        cdef double fitness = 0
        cdef int node
        target = {n: [] for n in self.target}
        for index in range(self.target_count):
            if not self.solve(angles[:, index]):
                # Punishment
                return HUGE_VAL
            for node in self.target:
                target[node].append(self.data_dict[node, -1])
        cdef Coordinate[:] path1, path2
        for node in self.target:
            path1 = array(target[node], dtype=object)
            if self.shape_only or self.wavelet_mode:
                _normalization(path1, 1)
                if self.wavelet_mode:
                    fitness += _cmp_wavelet(_wavelet(path1), self.target[node])
                    continue
            path2 = self.target[node]
            for index in range(self.target_count):
                fitness += (<Coordinate>path1[index]).distance(path2[index])
        return fitness

    cpdef object result(self, double[:] v):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        cdef int index = 0
        cdef int a_index = 0
        cdef VPoint vpoint
        for m in self.mapping_list:
            if type(m) is int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[index], v[index + 1])
                index += 2
            elif type(m) is str and m.startswith('A'):
                self.polar_angles[a_index] = v[index]
                a_index += 1
                index += 1
            else:
                self.mapping[m] = v[index]
                index += 1
        self.solve(v[self.l_base:self.l_base + self.target_count])
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
