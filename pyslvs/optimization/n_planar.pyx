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
from libc.math cimport M_PI, INFINITY as INF, hypot, fabs, atan2, cos, sin, sqrt
from numpy import array, zeros, float64 as f64
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
    cdef double[:, :] dist
    with gil:
        dist = zeros((len(p1), 2), dtype=f64)
    dist[:] = p1[:]
    cdef double neg_dist = 0
    cdef double pos_dist = 0
    cdef int i
    cdef double val, x, y
    for i in range(len(p1)):
        val = p1[i, 1] - a * p1[i, 0] - b
        x = p1[i, 0] - x_mean
        y = p1[i, 1] - y_mean
        if val < 0:
            neg_dist += sqrt(x * x + y * y)
        elif val > 0:
            pos_dist += sqrt(x * x + y * y)
    x = fabs(v2[0])
    y = fabs(v2[1])
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
    cdef double x, y
    for i in range(len(p1)):
        x = (p1[i, 0] - x_mean)
        y = (p1[i, 1] - y_mean)
        cxx += x * x
        cyy += y * y
        cxy += x * y
    cxx /= len(p1)
    cyy /= len(p1)
    cxy /= len(p1)
    cdef double[:, :] v
    with gil:
        _, v = eig(array([[cxx, cxy], [cxy, cyy]]))
    # Calculate the orientation of the axes
    cdef double a, b
    x, y = axes_v(v[:, 1], v[:, 0], p1, x_mean, y_mean)
    a, b = axes_v(v[:, 0], v[:, 1], p1, x_mean, y_mean)
    cdef double[:, :] axes
    with gil:
        axes = array([[x, y], [a, b]])
    cdef double theta_1 = atan2(axes[0][1], axes[0][0])
    cdef double theta_2 = atan2(axes[1][1], axes[1][0])
    # Calculate the rotation matrix
    if theta_1 * theta_2 > 0:
        return min(theta_1, theta_2)
    elif theta_1 * theta_2 < 0:
        if fabs(theta_1) < M_PI / 2:
            return min(theta_1, theta_2)
        else:
            return max(theta_1, theta_2)
    else:
        return -M_PI / 2 if -M_PI / 2 in {theta_1, theta_2} else 0


def norm_pca(path):
    """Normalization function by PCA."""
    cdef double[:, :] path_m = array(path, dtype=f64)
    _norm_pca(path_m)
    return [(x, y) for x, y in path_m]


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
    cdef double c = cos(alpha)
    cdef double s = sin(alpha)
    # Normalized the path points
    for i in range(len(p1)):
        p1[i, 0] -= x_mean
        p1[i, 1] -= y_mean
    with gil:
        p1 = array(p1) @ array([[c, -s], [s, c]])
    cdef int ind = 0
    cdef double p1_max = -INF
    cdef double p1_min = INF
    cdef double d_min = INF
    cdef double d
    for i in range(len(p1)):
        d = hypot(p1[i, 0], p1[i, 1])
        if d < d_min:
            d_min = d
            ind = i
        if p1[i, 0] > p1_max:
            p1_max = p1[i, 0]
        if p1[i, 0] < p1_min:
            p1_min = p1[i, 0]
    for i in range(len(p1)):
        p1[i, 0] /= p1_max - p1_min
        p1[i, 1] /= p1_max - p1_min
    # Swap to the starting point of the path
    roll(p1, ind)


cdef double trapezoidal_camp(double[:] a, double[:] b):
    """Error comparison by trapezoidal rule."""
    cdef double area = 0
    cdef int i
    for i in range(1, len(a)):
        area += abs(a[i - 1] + a[i] - b[i - 1] + b[i]) / 2
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
        self.target = array(mech['target'])
        self.len = len(self.target)
        # TODO: Normalization
        self.lb = array([1e-5] * 4 + [0.])
        self.ub = array([5.] * 4 + [2. * M_PI])
        raise NotImplementedError

    cdef double fitness(self, double[:] v) nogil:
        """Generate linkage with 5 parameters."""
        cdef double[:, :] p = c_uniform_path(v[None, :], self.len)[0]
        # TODO: Normalization
        raise NotImplementedError

    cpdef object result(self, double[:] v):
        return "M[" + ", ".join([(<VPoint> vp).expr()
                                 for vp in uniform_expr(v)]) + "]"
