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
from numpy import (
    zeros, array, arange, interp, argmax, concatenate, float64 as np_float,
)
from numpy.core.multiarray import correlate
from libc.math cimport (
    cos, sin, atan2, isnan, log, INFINITY as INF, HUGE_VAL, M_PI,
)
from .expression cimport Coord, VJoint, VPoint, distance
from .metaheuristics.utility cimport ObjFunc
from .tinycadlib cimport radians, pxy, ppp, plap, pllp, plpp, palp
from .topo_config cimport (
    t_config, symbol_str, I_LABEL, A_LABEL, Expr, PXY, PPP, PLA, PLAP, PLLP,
    PLPP, PALP, EStack,
)
from .bfgs cimport SolverSystem

try:
    from scipy.signal import fftconvolve
except ImportError:
    fftconvolve = None


def norm_path(path, scale=1):
    """Python wrapper of normalization function."""
    cdef double[:, :] path_m = array(path, dtype=np_float)
    _norm(path_m, scale)
    return [(x, y) for x, y in path_m]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _norm(double[:, :] path, double scale):
    """Normalization implementation inplace."""
    cdef double[:] centre = zeros(2, dtype=np_float)
    cdef double x, y
    for x, y in path:
        centre[0] += x
        centre[1] += y
    cdef int end = len(path)
    centre[0] /= end
    centre[1] /= end
    cdef double[:] angle = zeros(end + 1, dtype=np_float)
    cdef double[:] length = zeros(end + 1, dtype=np_float)
    cdef int sp = 0
    cdef int i
    for i in range(len(path)):
        x = path[i, 0]
        y = path[i, 1]
        angle[i] = atan2(y - centre[1], x - centre[0])
        length[i] = distance(centre[0], centre[1], x, y)
        if length[i] > length[end]:
            length[end] = length[i]
            angle[end] = angle[i]
            sp = end
    _aligned(path, sp)
    cdef double[:] bound = array([INF, -INF], dtype=np_float)
    cdef double a
    for i in range(len(path)):
        a = angle[i] - angle[end]
        path[i, 0] = length[i] * cos(a)
        path[i, 1] = length[i] * sin(a)
        if path[i, 0] < bound[0]:
            bound[0] = path[i, 0]
        if path[i, 0] > bound[1]:
            bound[1] = path[i, 0]
    scale /= (bound[1] - bound[0])
    _mul1d(path[:, 0], scale)
    _mul1d(path[:, 1], scale)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _aligned(double[:, :] path, int sp):
    """Split 1D path from sp, concatenate to end."""
    if sp == 0:
        return
    cdef double[:, :] tmp = zeros((sp, 2), dtype=np_float)
    cdef int i
    for i in range(sp):
        tmp[i, 0] = path[i, 0]
        tmp[i, 1] = path[i, 1]
    for i in range(sp, len(path)):
        path[i - sp, 0] = path[i, 0]
        path[i - sp, 1] = path[i, 1]
    for i in range(sp):
        path[len(path) - sp + i, 0] = tmp[i, 0]
        path[len(path) - sp + i, 1] = tmp[i, 1]


def curvature(path):
    r"""Calculate the signed curvature and return as an array.

    $$
    \kappa(t) = \frac{x'y'' - x''y'}{(x'^2 + y'^2)^\frac{3}{2}}
    $$
    """
    cdef double[:, :] path_m = array(path, dtype=np_float)
    return array(_curvature(path_m))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] _curvature(double[:, :] path):
    """Calculate the signed curvature."""
    cdef double[:, :] p1d = _derivative(path)
    cdef double[:, :] p2d = _derivative(p1d)
    cdef double[:] k = zeros(len(path) - 2, dtype=np_float)
    cdef int i
    for i in range(len(path) - 2):
        k[i] = ((p1d[i, 0] * p2d[i, 1] - p2d[i, 0] * p1d[i, 1])
                / (p1d[i, 0] * p1d[i, 0] + p1d[i, 1] * p1d[i, 1]) ** 1.5)
    return k


def derivative(double[:, :] p):
    """Differential function. Return $p'$."""
    return array(_derivative(p))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _derivative(double[:, :] p):
    """Differential function backend."""
    cdef double[:, :] pd = zeros((len(p), 2), dtype=np_float)
    cdef double max0 = 0
    cdef double max1 = 0
    cdef int i, j
    for i in range(len(p)):
        j = i + 1
        if j >= len(p):
            j = 0
        pd[i, 0] = p[j, 0] - p[i, 0]
        pd[i, 1] = p[j, 1] - p[i, 1]
        if pd[i, 0] > max0:
            max0 = pd[i, 0]
        if pd[i, 1] > max1:
            max1 = pd[i, 1]
    i = len(p) - 1
    if max0 == pd[i, 0] or max1 == pd[i, 1]:
        return pd[:i]
    else:
        return pd


def path_signature(double[:] k, double maximum = 100):
    r"""Require a curvature, return path signature.
    It's composed by curvature $\kappa$ and a $K$ value.

    $$
    K = \int^t_0 |\kappa(t)| dt
    $$

    >>> path_signature(curvature(...))
    """
    return array(_path_signature(k, maximum))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _path_signature(double[:] k, double maximum):
    """Require a curvature, return path signature."""
    cdef double[:, :] s = zeros((len(k), 2), dtype=np_float)
    cdef double interval = maximum / len(k)
    cdef double v = 0
    cdef int i
    for i in range(len(k)):
        s[i] = v
        v += interval
    s[:, 1] = k
    return s


def cross_correlation(double[:, :] p1, double[:, :] p2, double t = 0.1):
    r"""Compare signature and return as an 1d array.

    $$
    \begin{aligned}
    C_n(j, W, P) &= \left|\sum_i^{l_P} \frac{(W_{i + j}
    - \overline{W}_{j\rightarrow j + l_P})(P_i-\overline{P})}{
    \sqrt{\sum_i^{l_P}(W_{i + j} - \overline{W}_{j\rightarrow j + l_P})^2
    \sum_i^{l_P}(P_i - \overline{P})^2}}\right|
    \\
    S &= \arg\max\{C_n(j)\} t
    \end{aligned}
    $$

    >>> ps1 = path_signature(curvature(...))
    >>> ps2 = path_signature(curvature(...))
    >>> cc = cross_correlation(ps1, ps2)
    """
    return array(_cross_correlation(p1, p2, t), dtype=np_float)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] _cross_correlation(double[:, :] ps1, double[:, :] ps2, double t):
    """Compare path signature."""
    cdef double[:] p1 = interp(arange(0, _extr1d(ps1[:, 0], 1), t),
                               ps1[:, 0], ps1[:, 1])
    cdef double[:] p2 = interp(arange(0, _extr1d(ps2[:, 0], 1), t),
                               ps2[:, 0], ps2[:, 1])
    p1 = concatenate((p1, p1))
    cdef int s1 = len(p1)
    cdef int s2 = len(p2)
    cdef int n = s1 + s2 - 1
    if s2 >= s1:
        s1, s2 = s2, s1
    if (fftconvolve is not None and 1.89095737e-9 * 3 * n * log(n)
        < 2.1364985e-10 * (s1 - s2 + 1) * s2 - 1e-3
    ):
        return fftconvolve(p1, p2[::-1], mode='valid')
    else:
        return correlate(p1, p2)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _extr1d(double[:] s, bint op) nogil:
    """Max value of 1D slice."""
    cdef double m = -INF if op else INF
    for v in s:
        if op:
            if v > m:
                m = v
        elif v < m:
            m = v
    return m


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _mul1d(double[:] s, double v) nogil:
    """Inplace multiplication assignment of 1D slice."""
    cdef int i
    for i in range(len(s)):
        s[i] *= v


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _slice_nan2d(double[:, :] s):
    """Slice continuous view without NaN."""
    cdef int first = -1
    cdef int second = -1
    cdef int i
    cdef double v
    for i in range(len(s)):
        v = s[i, 1]
        if i == first + 1 and isnan(v):
            # The first of slice
            first = i
        elif isnan(v):
            # The end of slice
            second = i
            break
    if first == -1:
        first = 0
    if second == -1:
        second = len(s) - 1
    return s[first:second]


@cython.final
@cython.boundscheck(False)
@cython.wraparound(False)
cdef class FMatch(ObjFunc):
    """This class is used to verified kinematics of the linkage mechanism.

    A fast matching method that adds mapping angles to variables.
    """
    cdef bint bfgs_mode, shape_only, use_curvature, full_path, ordered
    cdef int target_count, input_count, l_base
    cdef list vpoints, mapping_list
    cdef dict placement, target, mapping, mapping_r, data_dict
    cdef double[:] polar_angles
    cdef EStack exprs

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
        #     'use_curvature': bool,
        #     'full_path': bool,
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
        self.use_curvature = mech.get('use_curvature', False)
        self.full_path = mech.get('full_path', False)
        if self.use_curvature:
            self.shape_only = False
        cdef int i, j
        cdef double[:, :] path
        for i in target:
            path = array(target[i], dtype=np_float)
            for j in range(i):
                if j in same:
                    i -= 1
            if self.shape_only:
                _norm(path, 1)
            if self.use_curvature:
                path = _path_signature(_curvature(path), 100)
            self.target[i] = path
        # Expressions
        self.vpoints = list(mech.get('expression', []))
        inputs = OrderedDict(mech.get('input', {}))
        self.input_count = len(inputs)
        status = {}
        self.exprs = t_config(self.vpoints, tuple(inputs.keys()), status)
        self.bfgs_mode = not all(status.values())
        # Data mapping
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_r = {v: k for k, v in self.mapping.items()}
        self.mapping_list = []
        # Bounds
        ub = []
        lb = []
        if len(ub) != len(lb):
            raise ValueError("upper and lower should be in the same size")
        # Position
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            ub.append(x + r)
            ub.append(y + r)
            lb.append(x - r)
            lb.append(y - r)
            self.mapping_list.append(i)
        cdef int p_base = len(ub)
        # Length of links
        link_upper = float(mech.get('upper', 100))
        link_lower = float(mech.get('lower', 0))
        i = 0
        cdef Expr expr
        for expr in self.exprs.stack:
            ub.append(link_upper)
            lb.append(link_lower)
            sym = frozenset({expr.c1.second, expr.target.second})
            self.mapping[sym] = None
            self.mapping_list.append(sym)
            if expr.func in {PLA, PLAP}:
                if expr.v2.first == A_LABEL:
                    ub.append(2 * M_PI)
                    lb.append(0)
                    sym = f"A{i}"
                    self.mapping[sym] = None
                    self.mapping_list.append(sym)
                    i += 1
            elif expr.func == PLLP:
                ub.append(link_upper)
                lb.append(link_lower)
                sym = frozenset({expr.c2.second, expr.target.second})
                self.mapping[sym] = None
                self.mapping_list.append(sym)
        self.l_base = len(ub)
        if self.use_curvature and self.full_path:
            # Scale factor
            ub.append(1)
            lb.append(1e-12)
        else:
            # Input nodes
            for start, end in inputs.values():
                ub.append(radians(start))
                lb.append(radians(end))
            # Angle rage
            ub[self.l_base:] *= self.target_count
            lb[self.l_base:] *= self.target_count
        self.ub = array(ub, dtype=np_float)
        self.lb = array(lb, dtype=np_float)
        # Swap upper and lower bound if reversed
        for i in range(len(self.ub)):
            if self.ub[i] < self.lb[i]:
                self.ub[i], self.lb[i] = self.lb[i], self.ub[i]
        # Allocate memory
        self.polar_angles = zeros(self.l_base - p_base, dtype=np_float)
        # Result list
        self.data_dict = {}

    cpdef bint is_two_kernel(self):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        return self.bfgs_mode

    cdef inline double get_len(self, str expr1, str expr2):
        return self.mapping[frozenset({self.mapping_r[expr1], self.mapping_r[
            expr2]})]

    cdef bint solve(self, double[:] input_list):
        self.data_dict.clear()
        cdef int i
        cdef VPoint vpoint
        cdef Coord coord1, coord2, coord3
        for i, vpoint in enumerate(self.vpoints):
            coord1 = Coord.__new__(Coord, vpoint.c[0, 0], vpoint.c[0, 1])
            if vpoint.type == VJoint.R:
                self.data_dict[self.mapping[i]] = coord1
                self.data_dict[i, -1] = coord1
            else:
                coord2 = Coord.__new__(Coord, vpoint.c[1, 0], vpoint.c[1, 1])
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
            if expr.func == PXY:
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord3 = pxy(
                    coord1,
                    vpoint.c[0, 0] - coord1.x,
                    vpoint.c[0, 1] - coord1.y
                )
            elif expr.func == PPP:
                coord3 = ppp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.data_dict[symbol_str(expr.c2)],
                    self.data_dict[symbol_str(expr.c3)],
                )
            elif expr.func in {PLA, PLAP}:
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
            elif expr.func == PALP:
                angle = self.polar_angles[a]
                a += 1
                coord3 = palp(
                    self.data_dict[symbol_str(expr.c1)],
                    angle,
                    self.get_len(symbol_str(expr.c2), target),
                    self.data_dict[symbol_str(expr.c2)],
                    expr.op
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
        # Solve
        try:
            solved_bfgs = SolverSystem(self.vpoints, {}, p_data_dict).solve()
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
                self.data_dict[i, -1] = Coord.__new__(
                    Coord,
                    solved_bfgs[i][0],
                    solved_bfgs[i][1]
                )
            else:
                self.data_dict[i, -1] = Coord.__new__(
                    Coord,
                    solved_bfgs[i][0][0],
                    solved_bfgs[i][0][1]
                )
                self.data_dict[i, -2] = Coord.__new__(
                    Coord,
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
        # Copy data
        cdef int i = 0
        cdef int j = 0
        for m in self.mapping_list:
            if type(m) is int:
                (<VPoint>self.vpoints[m]).locate(v[i], v[i + 1])
                i += 2
            elif type(m) is str and m.startswith('A'):
                self.polar_angles[j] = v[i]
                j += 1
                i += 1
            else:
                self.mapping[m] = v[i]
                i += 1
        # Solver
        # Angles (node, path length)
        cdef double[:, :] angles
        if self.use_curvature:
            angles = zeros((self.input_count, 36), dtype=np_float)
            for i in range(self.input_count):
                for j in range(0, 360, 10):
                    angles[i, j] = radians(j)
        else:
            angles = zeros((self.input_count, self.target_count),
                           dtype=np_float)
            for i in range(self.input_count):
                j = i + self.l_base
                angles[i, :] = v[j:j + self.target_count]
        cdef double fitness = 0
        cdef int node
        cdef Coord c
        target = {n: [] for n in self.target}
        for i in range(self.target_count):
            if not self.solve(angles[:, i]):
                # Punishment
                return HUGE_VAL
            for node in self.target:
                c = self.data_dict[node, -1]
                target[node].append((c.x, c.y))
        # Compare
        cdef double scale
        cdef double[:, :] path1, path2
        for node in self.target:
            path1 = array(target[node], dtype=np_float)
            path2 = array(self.target[node])
            if self.use_curvature:
                path1 = _slice_nan2d(path1)
                if len(path1) == 0:
                    return HUGE_VAL
                path1 = _path_signature(_curvature(path1), 100)
                scale = 1 / v[len(v) - 1]
                if not self.full_path:
                    scale *= _extr1d(path1[:, 0], 1) / _extr1d(path2[:, 0], 1)
                _mul1d(path2[:, 0], scale)
                j = argmax(_cross_correlation(path2, path1, 0.1))
                for i in range(len(path2)):
                    path2[i, 0] += j
                for i in range(self.target_count):
                    fitness += path1[i, 0] - path2[i, 0]
            else:
                if self.shape_only:
                    _norm(path1, 1)
                for i in range(self.target_count):
                    fitness += distance(path1[i, 0], path1[i, 1],
                                        path2[i, 0], path2[i, 1])
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
        cdef Coord coord1, coord2
        for i in range(len(self.vpoints)):
            vpoint = self.vpoints[i]
            coord1 = self.data_dict[i, -1]
            vpoint.locate(coord1.x, coord1.y)
            if vpoint.type != VJoint.R:
                coord2 = self.data_dict[i, -2]
                vpoint.move((coord1.x, coord1.y), (coord2.x, coord2.y))
            expressions.append(vpoint.expr())
        return "M[" + ", ".join(expressions) + "]"
