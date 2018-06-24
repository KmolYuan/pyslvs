# -*- coding: utf-8 -*-
#cython: language_level=3
# distutils: include_dirs = ./bfgs_solver
# distutils: sources = constraint_func.cpp, derivatives.cpp, solve.cpp

"""Wrapper of BFGS algorithm."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

cdef extern from "solve.h":
    
    cdef int Rough
    cdef int Fine
    cdef int MaxIterations
    
    cdef int Succsess
    cdef int NoSolution
    
    struct Point:
        double *x
        double *y
    
    struct Line:
        Point *p1
        Point *p2
    
    struct Constraint:
        pass
    
    Constraint PointOnPointConstraint(Point *, Point *)
    Constraint P2PDistanceConstraint(Point *, Point *, double *)
    Constraint HorizontalConstraint(Line *)
    
    int solve(double **, int, Constraint *, int, int)
    void derivatives(double **, double *, int, Constraint *, int)

from libc.stdlib cimport malloc, free


cpdef test_kernel():
    cdef int i
    cdef list params = [
        0, 0,
        5, 0,
        6, 5,
        6, 5,
    ]
    cdef double *parameters = <double *>malloc(len(params) * sizeof(double))
    for i in range(len(params)):
        parameters[i] = params[i]
    cdef double **pparameters = <double **>malloc(len(params) * sizeof(double *))
    for i in range(8):
        pparameters[i] = &parameters[i]
    cdef double constants[3]
    constants = [30, 10, 24]
    
    cdef Point points[5]
    points = [
        [pparameters[0], pparameters[1]],
        [pparameters[2], pparameters[3]],
        [pparameters[4], pparameters[5]],
        [pparameters[6], pparameters[7]],
        [&constants[0], &constants[1]],
    ]
    cdef Line lines[2]
    lines = [
        [&points[0], &points[1]],
        [&points[1], &points[2]],
    ]
    cdef Constraint cons[4]
    cons = [
        HorizontalConstraint(&lines[0]),
        PointOnPointConstraint(&points[2], &points[4]),
        PointOnPointConstraint(&points[3], &points[4]),
        P2PDistanceConstraint(&points[1], &points[2], &constants[2]),
    ]
    
    for i in range(5):
        print("Point {}: ({}, {})".format(i, points[i].x[0], points[i].y[0]))
    
    if solve(pparameters, len(params), cons, len(cons), Rough) == Succsess:
        print("A good Solution was found.")
    else:
        print("No valid Solutions were found from this start point.")
    
    for i in range(5):
        print("Point {}: ({}, {})".format(i, points[i].x[0], points[i].y[0]))
    
    cdef list zeros = [0] * 5
    cdef double *gradF = <double *>malloc(len(zeros) * sizeof(double))
    for i in range(len(zeros)):
        gradF[i] = zeros[i]
    derivatives(pparameters, gradF, len(zeros), cons, len(cons))
    for i in range(len(zeros)):
        print("GradF[{}]: {}".format(i, gradF[i]))
    
    free(parameters)
    free(pparameters)
    free(gradF)
