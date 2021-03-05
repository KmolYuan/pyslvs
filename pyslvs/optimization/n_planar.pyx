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
from libc.math cimport (
    M_PI, INFINITY as INF, HUGE_VAL, isnan, hypot, atan2, cos, sin, sqrt,
)
from numpy import array, float64 as f64, hstack
from numpy.linalg import eig
from numpy.fft import fft
from pyslvs.expression cimport VPoint
from pyslvs.tinycadlib cimport c_uniform_path, uniform_expr
from pyslvs.metaheuristics.utility cimport ObjFunc
from .f_planar cimport roll


cdef (double, double) axes_v(double[:] v1, double[:] v2, double[:, :] p1,
                             double x_mean, double y_mean) nogil:
    """Calculate the orientation vector of the axes."""
    cdef double a = v1[1] / v1[0]
    cdef double b = y_mean - a * x_mean
    cdef double neg_x = 0
    cdef double neg_y = 0
    cdef double pos_x = 0
    cdef double pos_y = 0
    cdef int i
    cdef double val, x, y
    for i in range(len(p1)):
        val = p1[i, 1] - a * p1[i, 0] - b
        x = p1[i, 0] - x_mean
        y = p1[i, 1] - y_mean
        if val < 0:
            neg_x += x * x
            neg_y += y * y
        elif val > 0:
            pos_x += x * x
            pos_y += y * y
    cdef double neg_dist = sqrt(neg_x) + sqrt(neg_y)
    cdef double pos_dist = sqrt(pos_x) + sqrt(pos_y)
    x = abs(v2[0])
    y = abs(v2[1])
    if a < 0:
        if neg_dist > pos_dist:
            return -x, -y
        else:
            return x, y
    else:
        if neg_dist > pos_dist:
            return x, -y
        else:
            return -x, y


cdef double rotation_angle(double[:, :] p1, double x_mean, double y_mean) nogil:
    """Calculate rotation angle."""
    # Calculate the principal component axes (eigenvectors)
    cdef double cxx = 0
    cdef double cyy = 0
    cdef double cxy = 0
    cdef int i
    cdef double x, y
    for i in range(len(p1)):
        x = p1[i, 0] - x_mean
        y = p1[i, 1] - y_mean
        cxx += x * x
        cyy += y * y
        cxy += x * y
    cxx /= len(p1)
    cyy /= len(p1)
    cxy /= len(p1)
    cdef double[:, :] v
    with gil:
        _, v = eig(array([[cxx, cxy], [cxy, cyy]], dtype=f64))
    # Calculate the orientation of the axes
    x, y = axes_v(v[:, 1], v[:, 0], p1, x_mean, y_mean)
    cdef double a1 = atan2(y, x)
    x, y = axes_v(v[:, 0], v[:, 1], p1, x_mean, y_mean)
    cdef double a2 = atan2(y, x)
    # Calculate the rotation matrix
    if a1 * a2 > 0:
        return min(a1, a2)
    elif a1 * a2 < 0:
        if abs(a1) < M_PI / 2:
            return min(a1, a2)
        else:
            return max(a1, a2)
    else:
        return -M_PI / 2 if -M_PI / 2 in {a1, a2} else 0


cdef double[:, :] rotate(double[:, :] p1, double a) nogil:
    cdef double c = cos(a)
    cdef double s = sin(a)
    with gil:
        return array(p1) @ array([[c, -s], [s, c]], dtype=f64)


def norm_pca(path):
    """Normalization function by PCA."""
    cdef double[:, :] path_m = array(path, dtype=f64)
    _norm_pca(path_m)
    return array(path_m)


cdef void _norm_pca(double[:, :] p1) nogil:
    """Normalization implementation."""
    cdef double x_mean = 0
    cdef double y_mean = 0
    cdef int i
    for i in range(len(p1)):
        x_mean += p1[i, 0]
        y_mean += p1[i, 1]
    x_mean /= len(p1)
    y_mean /= len(p1)
    cdef double alpha = rotation_angle(p1, x_mean, y_mean)
    # Normalized the path points
    for i in range(len(p1)):
        p1[i, 0] -= x_mean
        p1[i, 1] -= y_mean
    p1[:] = rotate(p1, alpha)
    cdef int ind = 0
    cdef double p1_max = -INF
    cdef double p1_min = INF
    cdef double d_min = INF
    cdef double d
    for i in range(len(p1)):
        d = hypot(p1[i, 0], p1[i, 1])
        if d < d_min:
            d_min = d
        if p1[i, 0] > p1_max:
            p1_max = p1[i, 0]
            if p1[i, 1] > 0:
                ind = i
        if p1[i, 0] < p1_min:
            p1_min = p1[i, 0]
    d_min = p1_max - p1_min
    for i in range(len(p1)):
        p1[i, 0] /= d_min
        p1[i, 1] /= d_min
    # Swap to the starting point of the path
    roll(p1, ind)


cdef void transform(double[:, :] target) nogil:
    """Apply normalization and FFT to target."""
    _norm_pca(target)
    cdef double[:] real, imag
    with gil:
        c = array(target, dtype=f64)
        c = c[:, 0] + c[:, 1] * 1j
        v = fft(hstack([c, c, c]))[len(target):len(target) * 2]
        real = v.real
        imag = v.imag
    target[:, 0] = real
    target[:, 1] = imag


cdef double trapezoidal_camp(double[:] a, double[:] b) nogil:
    """Error comparison by trapezoidal rule."""
    cdef double area = 0
    cdef int i
    cdef double right, left, total
    for i in range(1, len(a)):
        right = a[i] - b[i]
        left = a[i - 1] - b[i - 1]
        total = abs(right) + abs(left)
        if right * left < 0:
            area += (right * right + left * left) / total * 0.5
        else:
            area += total * 0.5
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
        self.target = array(mech['target'], dtype=f64)
        transform(self.target)
        self.len = len(self.target)
        self.lb = array([1e-5] * 4 + [0.], dtype=f64)
        self.ub = array([5.] * 4 + [2. * M_PI], dtype=f64)

    cdef double fitness(self, double[:] v) nogil:
        """Generate linkage with 5 parameters."""
        cdef double[:, :] p = c_uniform_path(v[None, :], self.len)[0]
        # NAN check
        cdef int i
        for i in range(len(p)):
            if isnan(p[i, 0]) or isnan(p[i, 1]):
                return HUGE_VAL
        transform(p)
        return (trapezoidal_camp(self.target[:, 0], p[:, 0]) +
                trapezoidal_camp(self.target[:, 1], p[:, 1]))

    cpdef object path(self, double[:] v):
        return array(c_uniform_path(v[None, :], self.len)[0])

    cpdef object result(self, double[:] v):
        return "M[" + ", ".join([(<VPoint> vp).expr()
                                 for vp in uniform_expr(v)]) + "]"
