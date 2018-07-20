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
from libc.math cimport (
    M_PI,
    cos,
    sin,
)
from libcpp.map cimport map
from sketch_solve cimport *
from tinycadlib cimport VPoint


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
    
    + vpoints: Sequence[VPoint]
    + inputs: [(b0, d0, a0), (b1, d1, a1), ...]
    """
    cdef dict vlinks = {}
    
    #Pre-count number of parameters.
    cdef int params_count = 0
    cdef int constants_count = 0
    cdef int slider_p_count = 0
    cdef int slider_rp_count = 0
    #sliders = {p_num: base_num}
    cdef map[int, int] sliders
    
    cdef int a = 0
    cdef int b = 0
    cdef int i
    cdef VPoint vpoint
    for i, vpoint in enumerate(vpoints):
        if vpoint.grounded():
            constants_count += 2
        else:
            params_count += 2
        if vpoint.type in {1, 2}:
            if vpoint.grounded():
                a += 1
            else:
                b += 1
            params_count += 4
            sliders[i] = slider_p_count + slider_rp_count
            if vpoint.type == 1:
                slider_p_count += 1
            else:
                slider_rp_count += 1
    
    cdef double *parameters = <double *>malloc(params_count * sizeof(double))
    cdef double **pparameters = <double **>malloc(params_count * sizeof(double *))
    cdef double *constants = <double *>malloc(constants_count * sizeof(double))
    cdef Point *points = <Point *>malloc(len(vpoints) * sizeof(Point))
    #Slider data
    cdef int slider_count = slider_p_count + slider_rp_count
    cdef Point *slider_bases = NULL
    cdef Point *slider_slots = NULL
    cdef Line *slider_lines = NULL
    cdef double *cons_angles = NULL
    if not sliders.empty():
        slider_bases = <Point *>malloc(slider_count * sizeof(Point))
        slider_slots = <Point *>malloc(slider_count * sizeof(Point))
    
    #Create parameters and link data.
    slider_count = 0
    a = 0
    b = 0
    cdef str vlink
    for i, vpoint in enumerate(vpoints):
        for vlink in vpoint.links:
            if vlink == 'ground':
                continue
            if vlink in vlinks:
                vlinks[vlink].append(i)
            else:
                vlinks[vlink] = [i]
        if vpoint.grounded():
            constants[b], constants[b + 1] = vpoint.c[0]
            if vpoint.type in {1, 2}:
                #Base point (slot) is fixed.
                slider_bases[slider_count] = [constants + b, constants + b + 1]
                #Slot point (slot) is moveable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a], pparameters[a + 1] = &parameters[a], &parameters[a + 1]
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
                #Pin is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
                pparameters[a], pparameters[a + 1] = &parameters[a], &parameters[a + 1]
                points[i] = [pparameters[a], pparameters[a + 1]]
                a += 2
            else:
                #Point is fixed.
                points[i] = [constants + b, constants + b + 1]
            b += 2
        else:
            if vpoint.type in {1, 2}:
                #Base point (slot) is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
                pparameters[a], pparameters[a + 1] = &parameters[a], &parameters[a + 1]
                slider_bases[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                #Slot point (slot) is moveable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a], pparameters[a + 1] = &parameters[a], &parameters[a + 1]
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
            #Point / Pin is moveable.
            parameters[a], parameters[a + 1] = vpoint.c[0]
            pparameters[a], pparameters[a + 1] = &parameters[a], &parameters[a + 1]
            points[i] = [pparameters[a], pparameters[a + 1]]
            a += 2
    
    cdef int cons_count = 0
    
    #Pre-count number of distence constraints.
    cdef int link_count
    for vlink in vlinks:
        link_count = len(vlinks[vlink])
        if link_count < 2:
            continue
        cons_count += 1 + 2 * (link_count - 2)
    cdef double *distences = <double *>malloc(cons_count * sizeof(double))
    
    #Pre-count number of slider constraints.
    slider_count = 0
    link_count = 0
    for a, b in sliders:
        link_count += 1
        cons_count += 2
        slider_count += 1
        if not vpoints[a].grounded():
            link_count += 1
        if vpoints[a].type == 1:
            cons_count += 1
            link_count += 1
            slider_count += 1
    if not sliders.empty():
        slider_lines = <Line *>malloc(link_count * sizeof(Line))
        cons_angles = <double *>malloc(slider_count * sizeof(double))
    
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
        cons[i] = P2PDistanceConstraint(
            slider_bases + sliders[a] if vpoints[a].is_slot_link(vlink) else points + a,
            slider_bases + sliders[b] if vpoints[b].is_slot_link(vlink) else points + b,
            distences + i
        )
        i += 1
        for c in vlinks[vlink][2:]:
            for d in (a, b):
                distences[i] = vpoints[c].distance(vpoints[d])
                cons[i] = P2PDistanceConstraint(
                    slider_bases + sliders[c] if vpoints[c].is_slot_link(vlink) else points + c,
                    slider_bases + sliders[d] if vpoints[d].is_slot_link(vlink) else points + d,
                    distences + i
                )
                i += 1
    
    #Add slider constraints.
    #i: Constraint number.
    #c: Slider line number.
    #d: Angle digit number.
    c = 0
    d = 0
    cdef int f1
    for a, b in sliders:
        slider_lines[c] = [slider_bases + b, slider_slots + b]
        if vpoints[a].grounded():
            #Slot grounded.
            cons_angles[d] = vpoints[a].angle / 180 * M_PI
            cons[i] = LineInternalAngleConstraint(slider_lines + c, cons_angles + d)
            i += 1
            cons[i] = PointOnLineConstraint(points + a, slider_lines + c)
        else:
            #Slider between links.
            try:
                vlink = vpoints[a].links[0]
                f1 = vlinks[vlink][0]
                if f1 == a:
                    f1 = vlinks[vlink][1]
            except IndexError:
                pass
            else:
                c += 1
                if vpoints[f1].is_slot_link(vlink):
                    slider_lines[c] = [points + a, slider_bases + sliders[f1]]
                else:
                    slider_lines[c] = [points + a, points + f1]
                cons_angles[d] = (vpoints[a].angle - vpoints[a].slope_angle(vpoints[f1])) / 180 * M_PI
                cons[i] = InternalAngleConstraint(
                    slider_lines + c - 1,
                    slider_lines + c,
                    cons_angles + d
                )
                i += 1
                cons[i] = PointOnLineConstraint(points + a, slider_lines + c - 1)
        d += 1
        c += 1
        i += 1
        #P joint.
        if vpoints[a].type == 1:
            try:
                vlink = vpoints[a].links[1]
                f1 = vlinks[vlink][0]
                if f1 == a:
                    f1 = vlinks[vlink][1]
            except IndexError:
                pass
            else:
                if vpoints[f1].is_slot_link(vlink):
                    slider_lines[c] = [points + a, slider_bases + sliders[f1]]
                else:
                    slider_lines[c] = [points + a, points + f1]
                cons_angles[d] = vpoints[a].slope_angle(vpoints[f1]) / 180 * M_PI - cons_angles[d - 1]
                cons[i] = InternalAngleConstraint(
                    slider_lines + c - (1 if vpoints[a].grounded() else 2),
                    slider_lines + c,
                    cons_angles + d
                )
                d += 1
                c += 1
                i += 1
    
    #Add angle constraints for input angles.
    #c: Input data count.
    c = 0
    cdef double angle
    for b, d, angle in inputs:
        lines[c] = [points + b, points + d]
        angles[c] = angle / 180 * M_PI
        cons[i] = LineInternalAngleConstraint(lines + c, angles + c)
        i += 1
        c += 1
    
    #Solve
    if solve(pparameters, params_count, cons, cons_count, Rough) != Succsess:
        raise Exception("No valid Solutions were found from this start point.")
    
    """Format:
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
            solved_points.append((
                (slider_bases[sliders[i]].x[0], slider_bases[sliders[i]].y[0]),
                (points[i].x[0], points[i].y[0])
            ))
    free(parameters)
    free(pparameters)
    free(constants)
    free(points)
    if not sliders.empty():
        free(slider_bases)
        free(slider_slots)
        free(slider_lines)
        free(cons_angles)
    free(distences)
    free(angles)
    free(lines)
    free(cons)
    return solved_points


cpdef void partial_solving(object vpoints, dict data_dict):
    """Partial constants solving.
    
    + vpoints: Sequence[VPoint]
    + TODO: Known coordinates import from data_dict.
    + data_dict: Dict[int, Tuple[float, float]]
    + targets: Set[int]
    """
    cdef int i
    cdef set targets = {i for i in range(len(vpoints)) if (i not in data_dict)}
    cdef set friends = set()
    cdef dict vlinks = {}
    
    #Pre-count number of parameters.
    cdef int point_count = 0
    cdef int params_count = 0
    cdef int constants_count = 0
    cdef int cons_count = 0
    cdef int slider_p_count = 0
    cdef int slider_rp_count = 0
    #sliders = {p_num: base_num}
    cdef map[int, int] sliders
    
    #Create vlink data.
    cdef int a = 0
    cdef int b = 0
    cdef VPoint vpoint
    for i, vpoint in enumerate(vpoints):
        for vlink in vpoint.links:
            if vlink == 'ground':
                continue
            if vlink in vlinks:
                vlinks[vlink].append(i)
            else:
                vlinks[vlink] = [i]
        if (i in targets) and (vpoint.type in {1, 2}):
            if vpoint.grounded():
                a += 1
            else:
                b += 1
            params_count += 4
            sliders[i] = slider_p_count + slider_rp_count
            if vpoint.type == 1:
                slider_p_count += 1
            else:
                slider_rp_count += 1
    
    #Counting params and constants.
    for vlink in tuple(vlinks):
        if not any([(i in targets) for i in vlinks[vlink]]):
            #Remove the link not contain any target.
            del vlinks[vlink]
            continue
        for i in vlinks[vlink]:
            point_count += 1
            if i in targets:
                params_count += 2
            else:
                friends.add(i)
                constants_count += 2
    
    cdef double *parameters = <double *>malloc(params_count * sizeof(double))
    cdef double **pparameters = <double **>malloc(params_count * sizeof(double *))
    cdef double *constants = <double *>malloc(constants_count * sizeof(double))
    cdef Point *points = <Point *>malloc(point_count * sizeof(Point))
    #Slider data
    cdef int slider_count = slider_p_count + slider_rp_count
    cdef Point *slider_bases = NULL
    cdef Point *slider_slots = NULL
    cdef Line *slider_lines = NULL
    cdef double *cons_angles = NULL
    if not sliders.empty():
        slider_bases = <Point *>malloc(slider_count * sizeof(Point))
        slider_slots = <Point *>malloc(slider_count * sizeof(Point))
