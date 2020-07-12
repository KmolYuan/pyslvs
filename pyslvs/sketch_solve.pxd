# -*- coding: utf-8 -*-
# cython: language_level=3

"""Wrapper header of BFGS algorithm.
If you need to use container to carry floating point data,
std::list is recommended instead of std::vector.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cdef extern from "bfgs_solver/solve.h" nogil:

    int Rough
    int Fine
    int MaxIterations

    int Success
    int NoSolution

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

    Constraint PointOnPointConstraint(Point *, Point *);
    Constraint P2PDistanceConstraint(Point *, Point *, double *)
    Constraint PointOnLineConstraint(Point *, Line *)
    Constraint InternalAngleConstraint(Line *, Line *, double *)
    Constraint LineInternalAngleConstraint(Line *, double *)

    int solve(double **, size_t, Constraint *, size_t, int)
