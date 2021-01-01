# -*- coding: utf-8 -*-
# cython: language_level=3

"""Wrapper header of BFGS algorithm.
If you need to use container to carry floating point data,
std::list is recommended instead of std::vector.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cdef extern from "bfgs_solver/solve.h" nogil:
    struct Point:
        double *x
        double *y

    struct Line:
        Point *p1
        Point *p2

    struct Constraint:
        unsigned type
        Point *point1
        Point *point2
        Line *line1
        Line *line2
        double *parameter

    Constraint point_on_point(Point *, Point *)
    Constraint p2p_distance(Point *, Point *, double *)
    Constraint point_on_line(Point *, Line *)
    Constraint internal_angle(Line *, Line *, double *)
    Constraint line_internal_angle(Line *, double *)
    bint solve(double **, size_t, Constraint *, size_t, bint)
