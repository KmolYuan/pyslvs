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


cpdef test_kernel():
    cdef double parameters[8]
    parameters = [
        0, 0,
        5, 0,
        6, 5,
        6, 5,
    ]
    cdef double *pparameters[8]
    cdef int i
    for i in range(8):
        pparameters[i] = &parameters[i]
    cdef double constants[3]
    constants = [30, 10, 24]
    
    cdef Point points[5]
    points[0].x = pparameters[0]
    points[0].y = pparameters[1]
    
    points[1].x = pparameters[2]
    points[1].y = pparameters[3]
    
    points[2].x = pparameters[4]
    points[2].y = pparameters[5]
    
    points[3].x = pparameters[6]
    points[3].y = pparameters[7]
    
    points[4].x = &constants[0]
    points[4].y = &constants[1]
    
    cdef Line lines[2]
    lines[0].p1 = &points[0]
    lines[0].p2 = &points[1]
    
    lines[1].p1 = &points[1]
    lines[1].p2 = &points[2]
    
    cdef Constraint cons[4]
    cons = [
        HorizontalConstraint(&lines[0]),
        PointOnPointConstraint(&points[2], &points[4]),
        PointOnPointConstraint(&points[3], &points[4]),
        P2PDistanceConstraint(&points[1], &points[2], &constants[2]),
    ]
    
    for i in range(5):
        print("Point {}: ({}, {})".format(i, points[i].x[0], points[i].y[0]))
    
    if solve(pparameters, 8, cons, 4, Rough) == Succsess:
        print("A good Solution was found.")
    else:
        print("No valid Solutions were found from this start point.")
    
    for i in range(5):
        print("Point {}: ({}, {})".format(i, points[i].x[0], points[i].y[0]))
    
    cdef double gradF[5]
    gradF = [0] * 5
    derivatives(pparameters, gradF, 5, cons, 4)
    for i in range(5):
        print("GradF[{}]: {}".format(i, gradF[i]))
