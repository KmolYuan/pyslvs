# -*- coding: utf-8 -*-
# cython: language_level=3
# distutils: include_dirs = ./bfgs_solver
# distutils: sources = constraint_func.cpp, derivatives.cpp, solve.cpp

"""Wrapper of BFGS algorithm."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from libc.stdlib cimport malloc, free
from tinycadlib cimport VPoint


cpdef void test_kernel():
    """Test kernel of Sketch Solve. Same as C++ version."""
    cdef int i
    cdef double constants[3]
    constants = [30, 10, 24]
    cdef list params = [
        0, 0,
        5, 0,
        6, 5,
        6, 5,
    ]
    cdef double *parameters = <double *>malloc(len(params) * sizeof(double))
    cdef double **pparameters = <double **>malloc(len(params) * sizeof(double *))
    for i in range(len(params)):
        #Just copy values, these are not their pointer.
        parameters[i] = params[i]
        #Pointer of parameter pointer.
        pparameters[i] = &parameters[i]
    
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
    
    cdef Point point
    for point in points:
        print("Point {}: ({}, {})".format(i, point.x[0], point.y[0]))
    
    if solve(pparameters, len(params), cons, len(cons), Rough) == Succsess:
        print("A good Solution was found.")
    else:
        print("No valid Solutions were found from this start point.")
    
    for point in points:
        print("Point {}: ({}, {})".format(i, point.x[0], point.y[0]))
    
    cdef int zeros = 5
    cdef double *gradF = <double *>malloc(zeros * sizeof(double))
    for i in range(zeros):
        gradF[i] = 0
    derivatives(pparameters, gradF, zeros, cons, len(cons))
    for i in range(zeros):
        print("GradF[{}]: {}".format(i, gradF[i]))
    
    free(parameters)
    free(pparameters)
    free(gradF)


cpdef list vpoint_solving(object vpoints, object angles):
    """Solving function from vpoint list.
    
    + angles: [(p0, link0_0, link0_1, a0), (p1, link1_0, link1_1, a1), ...]
    """
    cdef list params = []
    cdef list constants = []
    cdef list points = []
    cdef VPoint vpoint
    for vpoint in vpoints:
        if vpoint.grounded():
            constants.extend(vpoint.c[0])
        else:
            params.extend(vpoint.c[0])
    cdef double *parameters = <double *>malloc(len(params) * sizeof(double))
    cdef double **pparameters = <double **>malloc(len(params) * sizeof(double *))
    cdef int i
    for i in range(len(params)):
        #Just copy values, these are not their pointer.
        parameters[i] = params[i]
        #Pointer of parameter pointer.
        pparameters[i] = &parameters[i]
