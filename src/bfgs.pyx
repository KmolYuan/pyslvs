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
import numpy as np


cpdef tuple test_kernel():
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
    
    cdef int point_count = 5
    cdef list input_data = []
    cdef list output_data = []
    for i in range(point_count):
        input_data.append((points[i].x[0], points[i].y[0]))
    
    if solve(pparameters, len(params), cons, len(cons), Rough) != Succsess:
        raise Exception("No valid Solutions were found from this start point.")
    
    for i in range(point_count):
        output_data.append((points[i].x[0], points[i].y[0]))
    
    cdef double *gradF = <double *>malloc(point_count * sizeof(double))
    for i in range(point_count):
        gradF[i] = 0
    derivatives(pparameters, gradF, point_count, cons, len(cons))
    cdef list grad_data = []
    for i in range(point_count):
        grad_data.append(gradF[i])
    
    free(parameters)
    free(pparameters)
    free(gradF)
    return input_data, output_data, grad_data


cpdef list vpoint_solving(object vpoints, object inputs = []):
    """Solving function from vpoint list.
    
    + inputs: [(b0, d0, a0), (b1, d1, a1), ...]
    """
    cdef dict vlinks = {}
    
    #Pre-count number of parameters.
    cdef int params_count = 0
    cdef int constants_count = 0
    
    cdef VPoint vpoint
    for vpoint in vpoints:
        if vpoint.grounded():
            constants_count += 2
        else:
            params_count += 2
    
    cdef double *parameters = <double *>malloc(params_count * sizeof(double))
    cdef double **pparameters = <double **>malloc(params_count * sizeof(double *))
    cdef double *constants = <double *>malloc(constants_count * sizeof(double))
    cdef Point *points = <Point *>malloc(len(vpoints) * sizeof(Point))
    
    cdef int i
    cdef int a = 0
    cdef int b = 0
    cdef str vlink
    cdef Point point
    for i, vpoint in enumerate(vpoints):
        for vlink in vpoint.links:
            if vlink == 'ground':
                continue
            if vlink in vlinks:
                vlinks[vlink].append(i)
            else:
                vlinks[vlink] = [i]
        if not vpoint.grounded():
            parameters[a], parameters[a + 1] = vpoint.c[0]
            pparameters[a] = &parameters[a]
            pparameters[a + 1] = &parameters[a + 1]
            point = [pparameters[a], pparameters[a + 1]]
            a += 2
        else:
            constants[b], constants[b + 1] = vpoint.c[0]
            point = [constants + b, constants + b + 1]
            b += 2
        points[i] = point
    
    cdef int cons_count = 0
    
    #Pre-count number of distence constraints.
    cdef int link_count
    for vlink in vlinks:
        link_count = len(vlinks[vlink])
        if link_count < 2:
            continue
        cons_count += 1 + 2 * (link_count - 2)
    cdef double *distences = <double *>malloc(cons_count * sizeof(double))
    
    #Pre-count number of angle constraints.
    cdef int input_count = len(inputs)
    cdef double *angles = <double *>malloc(input_count * sizeof(double))
    cdef Line *lines = <Line *>malloc(input_count * sizeof(Line))
    cons_count += input_count
    
    #Pre-count number of constraints.
    cdef Constraint *cons = <Constraint *>malloc(cons_count * sizeof(Constraint))
    
    #Create distence constraints of each link.
    i = 0
    cdef int c, d
    for vlink in vlinks:
        link_count = len(vlinks[vlink])
        if link_count < 2:
            continue
        a = vlinks[vlink][0]
        b = vlinks[vlink][1]
        distences[i] = vpoints[a].distance(vpoints[b])
        cons[i] = P2PDistanceConstraint(points + a, points + b, distences + i)
        i += 1
        for c in vlinks[vlink][2:]:
            for d in (a, b):
                distences[i] = vpoints[c].distance(vpoints[d])
                cons[i] = P2PDistanceConstraint(points + c, points + d, distences + i)
                i += 1
    
    #Add angle constraints for input angles.
    link_count = 0
    cdef double angle
    for b, d, angle in inputs:
        lines[link_count] = [points + b, points + d]
        angles[link_count] = np.deg2rad(angle)
        cons[i] = LineInternalAngleConstraint(lines + link_count, angles + link_count)
        i += 1
        link_count += 1
    
    #Solve
    if solve(pparameters, params_count, cons, cons_count, Rough) != Succsess:
        raise Exception("No valid Solutions were found from this start point.")
    
    """
    Format:
    (R joint)
    [p0]: (p0_x, p0_y)
    [p1]: (p1_x, p1_y)
    (P or RP joint)
    [p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))
    """
    cdef list solved_points = []
    for i in range(len(vpoints)):
        if vpoints[i].type == 0:
            solved_points.append((points[i].x[0], points[i].y[0]))
        else:
            solved_points.append((vpoints[i].c[0], (points[i].x[0], points[i].y[0])))
    free(parameters)
    free(pparameters)
    free(constants)
    free(points)
    free(distences)
    free(cons)
    return solved_points
