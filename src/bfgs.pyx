# -*- coding: utf-8 -*-
# cython: language_level=3

"""Wrapper of BFGS algorithm.

Note of Pointer:
+ In Cython, pointer is more convenient then array.
    Because we can not "new" them or using "const" decorator on size_t.
+ There is NO pointer's "get value" operator in Cython,
    please use "index" operator.
+ Pointers can be plus with C's Integer, but not Python's.
    So please copy or declare to C's Integer.
"""

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
from cpython cimport bool
from sketch_solve cimport (
    Rough,
    Succsess,
    Point,
    Line,
    Constraint,
    HorizontalConstraint,
    PointOnPointConstraint,
    P2PDistanceConstraint,
    InternalAngleConstraint,
    PointOnLineConstraint,
    LineInternalAngleConstraint,
    solve,
    derivatives,
)
from pmks cimport VPoint


cdef inline double radians(double degree):
    """Deg to rad."""
    return degree / 180 * M_PI


cdef inline tuple sort_pair(tuple pair):
    """A sorted pair."""
    cdef int a, b
    a, b = pair
    return pair if a < b else (b, a)


cpdef list vpoint_solving(
    object vpoints,
    object inputs = [],
    dict data_dict = {}
):
    """Solving function from vpoint list.
    
    + vpoints: Sequence[VPoint]
    + inputs: [(b0, d0, a0), (b1, d1, a1), ...]
    
    Known coordinates import from data_dict.
    + data_dict: {0: (10.0, 20.0), ..., (0, 2): 30.0, ...}
    """
    #Sort pairs in data_dict.
    cdef object k
    for k in data_dict:
        if type(k) == tuple:
            data_dict[sort_pair(k)] = data_dict.pop(k)
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
        if vpoint.no_link():
            constants_count += 2
            continue
        if vpoint.grounded():
            constants_count += 2
        else:
            if i in data_dict:
                #Known coordinates.
                constants_count += 2
            else:
                params_count += 2
        if vpoint.type in {1, 2}:
            if i in data_dict:
                #Known coordinates are not slider.
                continue
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
    
    #Avoid no parameters.
    if not params_count:
        return []
    
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
        #Base point and slot point to determine a slot line.
        slider_bases = <Point *>malloc(slider_count * sizeof(Point))
        slider_slots = <Point *>malloc(slider_count * sizeof(Point))
    
    #Create parameters and link data.
    slider_count = 0
    a = 0
    b = 0
    cdef str vlink
    for i, vpoint in enumerate(vpoints):
        if vpoint.no_link():
            constants[b], constants[b + 1] = vpoint.c[0]
            points[i] = [constants + b, constants + b + 1]
            b += 2
            continue
        for vlink in vpoint.links:
            if vlink == 'ground':
                continue
            if vlink in vlinks:
                vlinks[vlink].append(i)
            else:
                vlinks[vlink] = [i]
        if vpoint.grounded():
            if i in data_dict:
                #Known coordinates.
                constants[b], constants[b + 1] = data_dict[i]
                points[i] = [constants + b, constants + b + 1]
                b += 2
                continue
            constants[b], constants[b + 1] = vpoint.c[0]
            if vpoint.type in {1, 2}:
                #Base point (slot) is fixed.
                slider_bases[slider_count] = [constants + b, constants + b + 1]
                #Slot point (slot) is moveable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a], pparameters[a + 1] = (parameters + a), (parameters + a + 1)
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
                #Pin is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[1]
                pparameters[a], pparameters[a + 1] = (parameters + a), (parameters + a + 1)
                points[i] = [pparameters[a], pparameters[a + 1]]
                a += 2
            else:
                #Point is fixed.
                points[i] = [constants + b, constants + b + 1]
            b += 2
        else:
            if i in data_dict:
                #Known coordinates.
                constants[b], constants[b + 1] = data_dict[i]
                points[i] = [constants + b, constants + b + 1]
                b += 2
                continue
            if vpoint.type in {1, 2}:
                #Base point (slot) is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
                pparameters[a], pparameters[a + 1] = (parameters + a), (parameters + a + 1)
                slider_bases[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                #Slot point (slot) is moveable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a], pparameters[a + 1] = (parameters + a), (parameters + a + 1)
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
                #Pin is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[1]
            else:
                #Point is moveable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
            pparameters[a], pparameters[a + 1] = (parameters + a), (parameters + a + 1)
            points[i] = [pparameters[a], pparameters[a + 1]]
            a += 2
    
    #Pre-count number of distence constraints.
    cdef int cons_count = 0
    cdef int c, d
    for vlink in vlinks:
        if len(vlinks[vlink]) < 2:
            continue
        a = vlinks[vlink][0]
        b = vlinks[vlink][1]
        if (a not in data_dict) or (b not in data_dict):
            cons_count += 1
        for c in vlinks[vlink][2:]:
            if c in data_dict:
                #Known coordinates.
                continue
            cons_count += 2
    cdef double *distences = <double *>malloc(cons_count * sizeof(double))
    
    #Pre-count number of slider constraints.
    c = 0
    d = 0
    for a, b in sliders:
        c += 1
        cons_count += 2
        d += 1
        if not vpoints[a].grounded():
            c += 1
        if vpoints[a].type == 1:
            cons_count += 1
            c += 1
            d += 1
    if not sliders.empty():
        slider_lines = <Line *>malloc(c * sizeof(Line))
        cons_angles = <double *>malloc(d * sizeof(double))
    
    #Pre-count number of angle constraints.
    cdef int input_count = len(inputs)
    cdef double *angles = NULL
    cdef Line *lines = NULL
    if input_count:
        angles = <double *>malloc(input_count * sizeof(double))
        lines = <Line *>malloc(input_count * sizeof(Line))
        cons_count += input_count
    
    #Pre-count number of constraints.
    cdef Constraint *cons = <Constraint *>malloc(cons_count * sizeof(Constraint))
    
    #Create distence constraints of each link.
    i = 0
    cdef Point *p1
    cdef Point *p2
    cdef tuple pair
    for vlink in vlinks:
        if len(vlinks[vlink]) < 2:
            continue
        a = vlinks[vlink][0]
        b = vlinks[vlink][1]
        if (a not in data_dict) or (b not in data_dict):
            distences[i] = vpoints[a].distance(vpoints[b])
            if a in data_dict:
                p1 = points + a
            elif vpoints[a].is_slot_link(vlink):
                p1 = slider_bases + sliders[a]
            else:
                p1 = points + a
            if b in data_dict:
                p2 = points + b
            elif vpoints[b].is_slot_link(vlink):
                p2 = slider_bases + sliders[b]
            else:
                p2 = points + b
            cons[i] = P2PDistanceConstraint(p1, p2, distences + i)
            i += 1
        for c in vlinks[vlink][2:]:
            if c in data_dict:
                #Known coordinates.
                continue
            for d in (a, b):
                pair = sort_pair((c, d))
                if pair in data_dict:
                    distences[i] = data_dict[pair]
                else:
                    distences[i] = vpoints[c].distance(vpoints[d])
                if vpoints[c].is_slot_link(vlink):
                    p1 = slider_bases + sliders[c]
                else:
                    p1 = points + c
                if d in data_dict:
                    p2 = points + d
                elif vpoints[d].is_slot_link(vlink):
                    p2 = slider_bases + sliders[d]
                else:
                    p2 = points + d
                cons[i] = P2PDistanceConstraint(p1, p2, distences + i)
                i += 1
    
    #Add slider constraints.
    #i: Constraint number.
    #c: Slider line number.
    #d: Angle digit number.
    c = 0
    d = 0
    cdef int f1
    cdef Line *slider_slot
    for a, b in sliders:
        #Base slot.
        slider_lines[c] = [slider_bases + b, slider_slots + b]
        slider_slot = slider_lines + c
        if vpoints[a].grounded():
            #Slot is grounded.
            cons_angles[d] = radians(vpoints[a].angle)
            cons[i] = LineInternalAngleConstraint(slider_slot, cons_angles + d)
            i += 1
            cons[i] = PointOnLineConstraint(points + a, slider_slot)
        else:
            #Slider between links.
            try:
                vlink = vpoints[a].links[0]
                #A base link friend.
                f1 = vlinks[vlink][0]
                if f1 == a:
                    f1 = vlinks[vlink][1]
            except IndexError:
                pass
            else:
                #FIXME: There still have simulation error.
                c += 1
                if vpoints[f1].is_slot_link(vlink):
                    #f1 is a slider, and it is be connected with slot link.
                    slider_lines[c] = [slider_bases + b, slider_bases + sliders[f1]]
                else:
                    #f1 is a R joint or it is not connected with slot link.
                    slider_lines[c] = [slider_bases + b, points + f1]
                cons_angles[d] = radians(vpoints[a].slope_angle(vpoints[f1]) - vpoints[a].angle)
                cons[i] = InternalAngleConstraint(
                    slider_slot,
                    slider_lines + c,
                    cons_angles + d
                )
                i += 1
                cons[i] = PointOnLineConstraint(points + a, slider_slot)
        d += 1
        c += 1
        i += 1
        #P joint.
        if vpoints[a].type == 1:
            for vlink in vpoints[a].links[1:]:
                #A base link friend.
                f1 = vlinks[vlink][0]
                if f1 == a:
                    if len(vlinks[vlink]) < 2:
                        #If no any friend.
                        continue
                    f1 = vlinks[vlink][1]
                if vpoints[f1].is_slot_link(vlink):
                    #f1 is a slider, and it is be connected with slot link.
                    slider_lines[c] = [points + a, slider_bases + sliders[f1]]
                else:
                    #f1 is a R joint or it is not connected with slot link.
                    slider_lines[c] = [points + a, points + f1]
                cons_angles[d] = radians(vpoints[a].slope_angle(vpoints[f1]) - vpoints[a].angle)
                cons[i] = InternalAngleConstraint(
                    slider_slot,
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
        angles[c] = radians(angle)
        cons[i] = LineInternalAngleConstraint(lines + c, angles + c)
        i += 1
        c += 1
    
    #Solve
    if solve(pparameters, params_count, cons, cons_count, Rough) != Succsess:
        raise Exception("No valid Solutions were found from this start point.")
    
    """solved_points: List[
        #R joint
        [p1]: (p1_x, p1_y)
        #P or RP joint
        [p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))
    ]
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
    free(slider_bases)
    free(slider_slots)
    free(slider_lines)
    free(cons_angles)
    free(distences)
    free(angles)
    free(lines)
    free(cons)
    return solved_points
