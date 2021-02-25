# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Fast planar linkage synthesis. (defect allowable)

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from collections import OrderedDict
from logging import getLogger
from numpy import (
    zeros, array, arange, interp, argmax, concatenate, float64 as f64,
)
from numpy.core.multiarray import correlate
from libc.math cimport (
    cos, sin, atan2, isnan, log, INFINITY as INF, HUGE_VAL, M_PI,
)
from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.set cimport set as cset
from pyslvs.expression cimport Coord, VPoint, distance
from pyslvs.metaheuristics.utility cimport ObjFunc
from pyslvs.tinycadlib cimport (
    quick_solve, I_LABEL, A_LABEL, P_LABEL, PLA, PLAP, PLLP, preprocessing,
    CCoord, Sym, SwappablePair,
)
from pyslvs.topo_config cimport t_config, EStack, Expr
from pyslvs.bfgs cimport SolverSystem

try:
    from scipy.signal import fftconvolve
except ImportError:
    fftconvolve = None

logger = getLogger()


def norm_path(path, scale=1):
    """Normalization function."""
    cdef double[:, :] path_m = array(path, dtype=f64)
    _norm(path_m, scale)
    return array(path_m)


cdef void _norm(double[:, :] path, double scale) nogil:
    """Normalization implementation."""
    cdef CCoord centre = CCoord(0, 0)
    cdef double x, y
    cdef int i
    for i in range(len(path)):
        centre.x += path[i, 0]
        centre.y += path[i, 1]
    cdef int end = len(path)
    centre.x /= end
    centre.y /= end
    cdef vector[double] angle = vector[double](end + 1)
    cdef vector[double] length = vector[double](end + 1)
    cdef int sp = 0
    for i in range(len(path)):
        x = path[i, 0]
        y = path[i, 1]
        angle[i] = atan2(y - centre.y, x - centre.x)
        length[i] = distance(centre.x, centre.y, x, y)
        if length[i] > length[end]:
            length[end] = length[i]
            angle[end] = angle[i]
            sp = end
    roll(path, sp)
    cdef CCoord bound = CCoord(INF, -INF)
    cdef double a
    for i in range(len(path)):
        a = angle[i] - angle[end]
        path[i, 0] = length[i] * cos(a)
        path[i, 1] = length[i] * sin(a)
        if path[i, 0] < bound.x:
            bound.x = path[i, 0]
        if path[i, 0] > bound.y:
            bound.y = path[i, 0]
    scale /= (bound.y - bound.x)
    _mul1d(path[:, 0], scale)
    _mul1d(path[:, 1], scale)


cdef void roll(double[:, :] path, int ind) nogil:
    """Split 1D path from sp, concatenate to end."""
    if ind == 0:
        return
    cdef double[:, :] tmp
    with gil:
        tmp = zeros((ind, 2), dtype=f64)
    tmp[:] = path[:ind]
    path[:len(path) - ind] = path[ind:]
    path[len(path) - ind:] = tmp[:]


def curvature(path):
    r"""Calculate the signed curvature and return as an array.

    $$
    \kappa(t) = \frac{x'y'' - x''y'}{(x'^2 + y'^2)^\frac{3}{2}}
    $$
    """
    cdef double[:, :] path_m = array(path, dtype=f64)
    return array(_curvature(path_m))


cdef double[:] _curvature(double[:, :] path) nogil:
    """Calculate the signed curvature."""
    cdef double[:, :] p1d = _derivative(path)
    cdef double[:, :] p2d = _derivative(p1d)
    cdef double[:] k
    with gil:
        k = zeros(len(path) - 2, dtype=f64)
    cdef int i
    for i in range(len(path) - 2):
        k[i] = ((p1d[i, 0] * p2d[i, 1] - p2d[i, 0] * p1d[i, 1])
                / (p1d[i, 0] * p1d[i, 0] + p1d[i, 1] * p1d[i, 1]) ** 1.5)
    return k


def derivative(double[:, :] p):
    """Differential function. Return $p'$."""
    return array(_derivative(p))


cdef double[:, :] _derivative(double[:, :] p) nogil:
    """Differential function backend."""
    cdef double[:, :] pd
    with gil:
        pd = zeros((len(p), 2), dtype=f64)
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

    >>> from pyslvs.optimization import curvature, path_signature
    >>> path_signature(curvature(...))
    """
    return array(_path_signature(k, maximum))


cdef double[:, :] _path_signature(double[:] k, double maximum) nogil:
    """Require a curvature, return path signature."""
    cdef double[:, :] s
    with gil:
        s = zeros((len(k), 2), dtype=f64)
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

    >>> from pyslvs.optimization import curvature, path_signature
    >>> ps1 = path_signature(curvature(...))
    >>> ps2 = path_signature(curvature(...))
    >>> from pyslvs.optimization import cross_correlation
    >>> cc = cross_correlation(ps1, ps2)
    """
    return array(_cross_correlation(p1, p2, t), dtype=f64)


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


cdef void _mul1d(double[:] s, double v) nogil:
    """Inplace multiplication assignment of 1D slice."""
    cdef int i
    for i in range(len(s)):
        s[i] *= v


cdef double[:, :] _slice_nan2d(double[:, :] s) nogil:
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
cdef class FPlanar(ObjFunc):
    """A fast matching method that adds mapping angles to variables.

    Allowing defects.
    """
    cdef bint bfgs_mode, shape_only, use_curvature, full_path, ordered
    cdef int target_count, target_len, input_count, l_base
    cdef public unsigned long callback
    cdef list vpoints
    cdef long[:] target_nodes, pivots
    cdef double[:, :, :] target
    cdef EStack exprs
    cdef cset[int] slider
    cdef map[Sym, CCoord] joint_pos
    cdef map[Sym, double] param

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
        self.target_len = check_set.pop()
        # Change the target paths into memory view
        self.target_count = len(target)
        self.target_nodes = zeros(self.target_count, dtype=int)
        self.target = zeros((self.target_count, self.target_len, 2), dtype=f64)
        same = mech.get('same', {})
        self.shape_only = mech.get('shape_only', False)
        self.use_curvature = mech.get('use_curvature', False)
        self.full_path = mech.get('full_path', False)
        if self.use_curvature:
            self.shape_only = False
        cdef int i, j
        cdef double[:, :] path
        for i, j in enumerate(target):
            self.target_nodes[i] = j
            path = array(target[j], dtype=f64)
            if self.shape_only:
                _norm(path, 1)
            if self.use_curvature:
                path = _path_signature(_curvature(path), 100)
            self.target[i, :, :] = path
        # Expressions (must be readonly)
        self.vpoints = list(mech.get('expression', []))
        self.pivots = array([i for i, vp in enumerate(self.vpoints)
                             if (<VPoint>vp).grounded()], dtype=int)
        self.slider = {i for i, vp in enumerate(self.vpoints)
                       if (<VPoint>vp).is_slider()}
        inputs = OrderedDict(mech.get('input', {}))
        self.input_count = len(inputs)
        status = {}
        self.exprs = t_config(self.vpoints, tuple(inputs.keys()), status)
        self.bfgs_mode = not all(status.values())
        if not preprocessing(self.exprs, self.vpoints, {p: 0. for p in inputs},
                             self.joint_pos, map[SwappablePair, double](),
                             self.param):
            raise ValueError("wrong number of input parameters")
        # Bounds
        ub = []
        lb = []
        # Position
        cdef double x, y, r
        for i in self.pivots:
            x, y, r = placement[i]
            ub.append(x + r)
            ub.append(y + r)
            lb.append(x - r)
            lb.append(y - r)
        # Length of links
        link_upper = float(mech.get('upper', 100))
        link_lower = float(mech.get('lower', 0))
        cdef Expr expr
        for expr in self.exprs.stack:
            ub.append(link_upper)
            lb.append(link_lower)
            if expr.func in {PLA, PLAP} and expr.v2.first == A_LABEL:
                # The included angle of the link
                ub.append(2 * M_PI)
                lb.append(0)
            elif expr.func == PLLP:
                ub.append(link_upper)
                lb.append(link_lower)
        # The start of the angle parameters
        self.l_base = len(ub)
        if self.use_curvature and self.full_path:
            # Scale factor
            ub.append(1)
            lb.append(1e-12)
        else:
            # Input nodes
            for start, end in inputs.values():
                ub.append(start / 180 * M_PI)
                lb.append(end / 180 * M_PI)
            # Angle rage (input count * target count)
            ub[self.l_base:] *= self.target_len
            lb[self.l_base:] *= self.target_len
        self.ub = array(ub, dtype=f64)
        self.lb = array(lb, dtype=f64)
        # Swap upper and lower bound if reversed
        for i in range(len(self.ub)):
            if self.ub[i] < self.lb[i]:
                self.ub[i], self.lb[i] = self.lb[i], self.ub[i]
        self.callback = 0

    cpdef bint is_two_kernel(self):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        return self.bfgs_mode

    cdef double fitness(self, double[:] v) nogil:
        """Return the difference of the path signature.

        + Position of fix joints.
        + Link lengths.
        + Angle corresponding to the target points.
        """
        self.callback += 1
        cdef map[int, vector[CCoord]] target
        cdef map[Sym, double] param = map[Sym, double](self.param)
        cdef map[Sym, CCoord] joint_pos = map[Sym, CCoord](self.joint_pos)
        cdef bint ok
        cdef int i, j, node, vi
        cdef double x, y
        cdef Expr expr
        cdef CCoord c
        cdef pair[Sym, CCoord] jp
        for i in range(self.target_len):
            vi = 0
            for node in self.pivots:
                joint_pos[Sym(P_LABEL, node)] = CCoord(v[vi], v[vi + 1])
                vi += 2
            # Input parameters (length)
            for expr in self.exprs.stack:
                param[expr.v1] = v[vi]
                vi += 1
                if expr.v2.first == I_LABEL:
                    continue
                param[expr.v2] = v[vi]
                vi += 1
            # Input parameters (angles)
            for vi in range(self.input_count):
                param[Sym(I_LABEL, vi)] = v[self.l_base + vi
                                            + i * self.target_count]
            # Solve
            ok, joint_pos = quick_solve(self.exprs.stack, joint_pos, param)
            if not ok:
                return HUGE_VAL
            if self.bfgs_mode:
                with gil:
                    data_dict = {}
                    for jp in joint_pos:
                        data_dict[jp.first.second] = Coord.__new__(Coord,
                                                                   jp.second.x,
                                                                   jp.second.y)
                    # Solve
                    try:
                        solved_bfgs = SolverSystem(self.vpoints, {}, data_dict).solve()
                    except ValueError:
                        return HUGE_VAL
            # Collecting
            for node in self.target_nodes:
                if self.bfgs_mode:
                    with gil:
                        x, y = solved_bfgs[node][0]
                    target[node].push_back(CCoord(x, y))
                else:
                    c = joint_pos[Sym(P_LABEL, node)]
                    target[node].push_back(c)
        # Compare
        cdef double fitness = 0
        cdef double scale
        cdef double[:, :] path1, path2
        for vi, node in enumerate(self.target_nodes):
            if self.use_curvature:
                with gil:
                    path1 = zeros((self.target_len, 2), dtype=f64)
                    path2 = array(self.target[vi], dtype=f64)
                for i in range(self.target_len):
                    c = target[node][i]
                    path1[i, 0] = c.x
                    path1[i, 1] = c.y
                path1 = _slice_nan2d(path1)
                if len(path1) == 0:
                    return HUGE_VAL
                path1 = _path_signature(_curvature(path1), 100)
                scale = 1 / v[len(v) - 1]
                if not self.full_path:
                    scale *= _extr1d(path1[:, 0], 1) / _extr1d(path2[:, 0], 1)
                _mul1d(path2[:, 0], scale)
                with gil:
                    j = argmax(_cross_correlation(path2, path1, 0.1))
                for i in range(len(path2)):
                    path2[i, 0] += j
                for i in range(self.target_len):
                    fitness += path1[i, 0] - path2[i, 0]
            else:
                if self.shape_only:
                    with gil:
                        path1 = zeros((self.target_len, 2), dtype=f64)
                    for i in range(self.target_len):
                        c = target[node][i]
                        path1[i, 0] = c.x
                        path1[i, 1] = c.y
                    _norm(path1, 1)
                    for i in range(self.target_len):
                        target[node][i].x = path1[i, 0]
                        target[node][i].y = path1[i, 1]
                for i in range(self.target_len):
                    c = target[node][i]
                    fitness += distance(c.x, c.y,
                                        self.target[vi, i, 0],
                                        self.target[vi, i, 1])
        return fitness

    cpdef object result(self, double[:] v):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        cdef map[Sym, double] param = map[Sym, double](self.param)
        cdef map[Sym, CCoord] joint_pos = map[Sym, CCoord](self.joint_pos)
        cdef int vi = 0
        cdef int node
        for node in self.pivots:
            joint_pos[Sym(P_LABEL, node)] = CCoord(v[vi], v[vi + 1])
            vi += 2
        cdef Expr expr
        for expr in self.exprs.stack:
            param[expr.v1] = v[vi]
            vi += 1
            if expr.v2.first == I_LABEL:
                continue
            param[expr.v2] = v[vi]
            vi += 1
        # Input parameters (angles)
        for vi in range(self.input_count):
            param[Sym(I_LABEL, vi)] = v[self.l_base + vi]
        # Solve
        _, joint_pos = quick_solve(self.exprs.stack, joint_pos, param)
        cdef pair[Sym, CCoord] jp
        if self.bfgs_mode:
            data_dict = {}
            for jp in joint_pos:
                data_dict[jp.first.second] = Coord.__new__(Coord,
                                                           jp.second.x,
                                                           jp.second.y)
            # Solve
            solved_bfgs = SolverSystem(self.vpoints, {}, data_dict).solve()
        # Collecting
        cdef CCoord c
        cdef VPoint vp
        cdef double x, y
        for node in range(len(self.vpoints)):
            vp = self.vpoints[node]
            if self.bfgs_mode:
                x, y = solved_bfgs[node][0]
                vp.locate(x, y)
                if vp.is_slider():
                    vp.move(solved_bfgs[node][0], solved_bfgs[node][1])
            else:
                c = joint_pos[Sym(P_LABEL, node)]
                if vp.is_slider():
                    vp.locate(vp.c[0, 0], vp.c[0, 1])
                    vp.move((vp.c[0, 0], vp.c[0, 1]), (c.x, c.y))
                else:
                    vp.locate(c.x, c.y)
        expr_str = "M[" + ", ".join([(<VPoint> vp).expr()
                                     for vp in self.vpoints]) + "]"
        logger.debug(f"Result: {expr_str}")
        return expr_str
