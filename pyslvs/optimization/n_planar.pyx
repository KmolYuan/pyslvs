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
from libc.math cimport M_PI, fabs, atan2, cos, sin
from numpy import array, mean, matmul, max, min, argmin, roll, sum, insert
from numpy.linalg import norm, eig
from numpy.fft import fft
from pyslvs.expression cimport VPoint
from pyslvs.tinycadlib cimport c_uniform_path, uniform_expr
from pyslvs.metaheuristics.utility cimport ObjFunc


def axes_v(v1, v2, p1, x_mean, y_mean):
    """Calculate the orientation vector of the axes."""
    a = v1[1] / v1[0]
    b = y_mean - a * x_mean
    val = p1[:, 1] - a * p1[:, 0] - b
    neg_dist = norm((p1[val < 0] - mean(p1, axis=0)), axis=1).sum()
    pos_dist = norm((p1[val > 0] - mean(p1, axis=0)), axis=1).sum()
    v2 = abs(v2)
    if a < 0:
        if neg_dist > pos_dist:
            axis = [-v2[0], -v2[1]]
        else:
            axis = [v2[0], v2[1]]
    else:
        if neg_dist > pos_dist:
            axis = [v2[0], -v2[1]]
        else:
            axis = [-v2[0], v2[1]]
    return axis


def rotation_angle(p1):
    """Calculate rotation angle."""
    # Calculate the mean and std of path points
    x_mean, y_mean = mean(p1, axis=0)
    # Calculate the principal component axes (eigenvectors)
    m = p1.shape[0]
    cxx = sum((p1[:, 0] - x_mean) ** 2) / m
    cyy = sum((p1[:, 1] - y_mean) ** 2) / m
    cxy = sum((p1[:, 0] - x_mean) * (p1[:, 1] - y_mean)) / m
    c = array(((cxx, cxy), (cxy, cyy)))
    w, v = eig(c)
    # Calculate the orientation of the axes
    axes = array([
        axes_v(v[:, 1], v[:, 0], p1, x_mean, y_mean),
        axes_v(v[:, 0], v[:, 1], p1, x_mean, y_mean),
    ])
    theta_1 = atan2(axes[0][1], axes[0][0])
    theta_2 = atan2(axes[1][1], axes[1][0])
    pca = insert(axes, 1, (0, 0), axis=0) + mean(p1, axis=0)
    # Calculate the rotation matrix
    if theta_1 * theta_2 > 0:
        alpha = theta_1 if theta_1 < theta_2 else theta_2
    elif theta_1 * theta_2 < 0:
        if fabs(theta_1) < M_PI / 2:
            alpha = theta_1 if theta_1 < theta_2 else theta_2
        else:
            alpha = theta_1 if theta_1 > theta_2 else theta_2
    else:
        alpha = -M_PI / 2 if theta_1 == -M_PI / 2 or theta_2 == -M_PI / 2 else 0
    return alpha, pca


def normalization(p1):
    """Normalization function."""
    alpha, pca = rotation_angle(p1)
    c = cos(alpha)
    s = sin(alpha)
    r = array([[c, s], [-s, c]])
    # Normalized the path points
    x_mean, y_mean = mean(p1, axis=0)
    p1_normal = matmul(r, (p1 - array([x_mean, y_mean])).T).T
    p1_normal /= max(p1_normal[:, 0]) - min(p1_normal[:, 0])
    # Find the starting point of the path
    ind = argmin(norm(p1_normal, axis=1))
    return roll(p1_normal, p1_normal.shape[0] - ind, axis=0)


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
