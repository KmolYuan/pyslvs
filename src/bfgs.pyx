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

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from libc.stdlib cimport malloc, free
from libc.math cimport (
    M_PI,
    cos,
    sin,
)
from libcpp.map cimport map
from sketch_solve cimport (
    Rough,
    Succsess,
    Point,
    Line,
    Constraint,
    PointOnPointConstraint,
    P2PDistanceConstraint,
    InternalAngleConstraint,
    PointOnLineConstraint,
    LineInternalAngleConstraint,
    solve,
)
from pmks cimport VJoint, VPoint


@cython.cdivision
cdef inline double _radians(double degree):
    """Deg to rad."""
    return degree / 180 * M_PI


cdef inline tuple _sort_pair(int a, int b):
    """A sorted pair."""
    return (a, b) if a < b else (b, a)


cpdef list vpoint_solving(
    object vpoints,
    object inputs = None,
    object data_dict = None
):
    """Solving function from vpoint list.

    + vpoints: Sequence[VPoint]
    + inputs: [(b0, d0, a0), (b1, d1, a1), ...]

    Known coordinates import from data_dict.
    + data_dict: {0: (10.0, 20.0), ..., (0, 2): 30.0, ...}
    """
    # Blank sequences.
    if inputs is None:
        inputs = []
    if data_dict is None:
        data_dict = {}

    # Sort pairs in data_dict.
    cdef object k
    for k in data_dict:
        if type(k) == tuple:
            data_dict[_sort_pair(k[0], k[1])] = data_dict.pop(k)
    cdef dict vlinks = {}

    # Pre-count number of parameters.
    cdef int params_count = 0
    cdef int constants_count = 0
    cdef int slider_p_count = 0
    cdef int slider_rp_count = 0
    # sliders = {p_num: base_num}
    cdef map[int, int] sliders

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
                # Known coordinates.
                constants_count += 2
            else:
                if vpoint.pin_grounded():
                    constants_count += 2
                else:
                    params_count += 2
        if vpoint.type in {VJoint.P, VJoint.RP}:
            if i in data_dict:
                # Known coordinates are not slider.
                continue
            params_count += 4
            sliders[i] = slider_p_count + slider_rp_count
            if vpoint.type == VJoint.P:
                slider_p_count += 1
            else:
                slider_rp_count += 1

    # Avoid no parameters.
    if not params_count:
        return []

    cdef double *parameters = <double *>malloc(params_count * sizeof(double))
    cdef double **pparameters = <double **>malloc(params_count * sizeof(double *))
    cdef double *constants = <double *>malloc(constants_count * sizeof(double))
    cdef Point *points = <Point *>malloc(len(vpoints) * sizeof(Point))
    # Slider data
    cdef int slider_count = slider_p_count + slider_rp_count
    cdef Point *slider_bases = NULL
    cdef Point *slider_slots = NULL
    if not sliders.empty():
        # Base point and slot point to determine a slot line.
        slider_bases = <Point *>malloc(slider_count * sizeof(Point))
        slider_slots = <Point *>malloc(slider_count * sizeof(Point))

    # Create parameters and link data.
    slider_count = 0
    cdef int a = 0
    cdef int b = 0
    cdef str vlink
    cdef double offset
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
                # Known coordinates.
                constants[b], constants[b + 1] = data_dict[i]
                points[i] = [constants + b, constants + b + 1]
                b += 2
                continue
            constants[b], constants[b + 1] = vpoint.c[0]
            if vpoint.type in {VJoint.P, VJoint.RP}:
                # Base point (slot) is fixed.
                slider_bases[slider_count] = [constants + b, constants + b + 1]
                # Slot point (slot) is movable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a] = parameters + a
                pparameters[a + 1] = parameters + a + 1
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
                # Pin is movable.
                parameters[a], parameters[a + 1] = vpoint.c[1]
                if vpoint.has_offset() and (vpoint.true_offset() <= 0.1):
                    offset = 0.1 if vpoint.offset() > 0 else -0.1
                    parameters[a] += offset
                    parameters[a + 1] += offset
                pparameters[a] = parameters + a
                pparameters[a + 1] = parameters + a + 1
                points[i] = [pparameters[a], pparameters[a + 1]]
                a += 2
            else:
                # Point is fixed.
                points[i] = [constants + b, constants + b + 1]
            b += 2
        else:
            if i in data_dict:
                # Known coordinates.
                constants[b], constants[b + 1] = data_dict[i]
                points[i] = [constants + b, constants + b + 1]
                b += 2
                continue
            if vpoint.type in {VJoint.P, VJoint.RP}:
                # Base point (slot) is movable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
                pparameters[a] = parameters + a
                pparameters[a + 1] = parameters + a + 1
                slider_bases[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                # Slot point (slot) is movable.
                parameters[a] = vpoint.c[0][0] + cos(vpoint.angle)
                parameters[a + 1] = vpoint.c[0][1] + sin(vpoint.angle)
                pparameters[a] = parameters + a
                pparameters[a + 1] = parameters + a + 1
                slider_slots[slider_count] = [pparameters[a], pparameters[a + 1]]
                a += 2
                slider_count += 1
                if vpoint.pin_grounded():
                    # Pin is fixed.
                    constants[b], constants[b + 1] = vpoint.c[1]
                    points[i] = [constants + b, constants + b + 1]
                else:
                    # Pin is movable.
                    parameters[a], parameters[a + 1] = vpoint.c[1]
                    if vpoint.has_offset() and (vpoint.true_offset() <= 0.1):
                        offset = 0.1 if vpoint.offset() > 0 else -0.1
                        parameters[a] += offset
                        parameters[a + 1] += offset
                    pparameters[a] = parameters + a
                    pparameters[a + 1] = parameters + a + 1
                    points[i] = [pparameters[a], pparameters[a + 1]]
                    a += 2
            else:
                # Point is movable.
                parameters[a], parameters[a + 1] = vpoint.c[0]
                pparameters[a] = parameters + a
                pparameters[a + 1] = parameters + a + 1
                points[i] = [pparameters[a], pparameters[a + 1]]
                a += 2

    # Pre-count number of distance constraints.
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
                # Known coordinates.
                continue
            cons_count += 2
    cdef double *distances = <double *>malloc(cons_count * sizeof(double))

    # Pre-count number of slider constraints.
    c = 0
    d = 0
    cdef int f1
    cdef int slider_offset_count = 0
    for a, b in sliders:
        c += 1
        d += 1
        if vpoints[a].grounded():
            cons_count += 2
            if vpoints[a].has_offset():
                cons_count += 1
                if vpoints[a].offset():
                    slider_offset_count += 1
        else:
            for vlink in vpoints[a].links[:1]:
                f1 = vlinks[vlink][0]
                if f1 == a:
                    if len(vlinks[vlink]) < 2:
                        # If no any friend.
                        continue
                c += 1
                cons_count += 2
                if vpoints[a].has_offset():
                    cons_count += 1
                    if vpoints[a].offset():
                        slider_offset_count += 1
        if vpoints[a].type == VJoint.P:
            cons_count += 1
            c += 1
            d += 1

    cdef Line *slider_lines = NULL
    cdef double *cons_angles = NULL
    cdef double *slider_offset = NULL
    if not sliders.empty():
        slider_lines = <Line *>malloc(c * sizeof(Line))
        cons_angles = <double *>malloc(d * sizeof(double))
    if slider_offset_count:
        slider_offset = <double *>malloc(slider_offset_count * sizeof(double))

    # Pre-count number of angle constraints.
    cdef int input_count = 0
    cdef double angle
    for b, d, angle in inputs:
        if b == d:
            continue
        input_count += 1

    cdef double *angles = NULL
    cdef Line *lines = NULL
    if input_count:
        angles = <double *>malloc(input_count * sizeof(double))
        lines = <Line *>malloc(input_count * sizeof(Line))
        cons_count += input_count

    # Pre-count number of constraints.
    cdef Constraint *cons = <Constraint *>malloc(cons_count * sizeof(Constraint))

    # Create distance constraints of each link.
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
            distances[i] = vpoints[a].distance(vpoints[b])
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
            cons[i] = P2PDistanceConstraint(p1, p2, distances + i)
            i += 1
        for c in vlinks[vlink][2:]:
            if c in data_dict:
                # Known coordinates.
                continue
            for d in (a, b):
                pair = _sort_pair(c, d)
                if pair in data_dict:
                    distances[i] = data_dict[pair]
                else:
                    distances[i] = vpoints[c].distance(vpoints[d])
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
                cons[i] = P2PDistanceConstraint(p1, p2, distances + i)
                i += 1

    # Add slider constraints.
    # i: Constraint number.
    # c: Slider line number.
    # d: Angle digit number.
    c = 0
    d = 0
    slider_offset_count = 0
    cdef Line *slider_slot
    for a, b in sliders:
        # Base point.
        vpoint = vpoints[a]
        p1 = points + a
        # Base slot.
        slider_lines[c] = [slider_bases + b, slider_slots + b]
        slider_slot = slider_lines + c
        if vpoint.grounded():
            # Slot is grounded.
            cons_angles[d] = _radians(vpoint.angle)
            cons[i] = LineInternalAngleConstraint(slider_slot, cons_angles + d)
            i += 1
            cons[i] = PointOnLineConstraint(p1, slider_slot)
            if vpoint.has_offset():
                i += 1
                if vpoint.offset():
                    slider_offset[slider_offset_count] = vpoint.offset()
                    cons[i] = P2PDistanceConstraint(
                        slider_bases + b,
                        p1,
                        slider_offset + slider_offset_count
                    )
                    slider_offset_count += 1
                else:
                    cons[i] = PointOnPointConstraint(slider_bases + b, p1)
        else:
            # Slider between links.
            for vlink in vpoint.links[:1]:
                # A base link friend.
                f1 = vlinks[vlink][0]
                if f1 == a:
                    if len(vlinks[vlink]) < 2:
                        # If no any friend.
                        continue
                    f1 = vlinks[vlink][1]
                c += 1
                if vpoints[f1].is_slot_link(vlink):
                    # f1 is a slider, and it is be connected with slot link.
                    p2 = slider_bases + sliders[f1]
                else:
                    # f1 is a R joint or it is not connected with slot link.
                    p2 = points + f1
                slider_lines[c] = [slider_bases + b, p2]
                cons_angles[d] = _radians(vpoint.slope_angle(vpoints[f1]) - vpoint.angle)
                cons[i] = InternalAngleConstraint(
                    slider_slot,
                    slider_lines + c,
                    cons_angles + d
                )
                i += 1
                cons[i] = PointOnLineConstraint(p1, slider_slot)

                if not vpoint.has_offset():
                    continue

                i += 1
                if vpoint.offset():
                    slider_offset[slider_offset_count] = vpoint.offset()
                    cons[i] = P2PDistanceConstraint(
                        slider_bases + b,
                        p1,
                        slider_offset + slider_offset_count
                    )
                    slider_offset_count += 1
                else:
                    cons[i] = PointOnPointConstraint(slider_bases + b, p1)
        d += 1
        c += 1
        i += 1
        if vpoint.type != VJoint.P:
            continue

        for vlink in vpoint.links[1:]:
            # A base link friend.
            f1 = vlinks[vlink][0]
            if f1 == a:
                if len(vlinks[vlink]) < 2:
                    # If no any friend.
                    continue
                f1 = vlinks[vlink][1]
            if vpoints[f1].is_slot_link(vlink):
                # f1 is a slider, and it is be connected with slot link.
                p2 = slider_bases + sliders[f1]
            else:
                # f1 is a R joint or it is not connected with slot link.
                p2 = points + f1
            slider_lines[c] = [p1, p2]
            cons_angles[d] = _radians(vpoint.slope_angle(vpoints[f1]) - vpoint.angle)
            cons[i] = InternalAngleConstraint(
                slider_slot,
                slider_lines + c,
                cons_angles + d
            )
            d += 1
            c += 1
            i += 1

    # Add angle constraints for input angles.
    # c: Input data count.
    c = 0
    for b, d, angle in inputs:
        if b == d:
            continue
        lines[c] = [points + b, points + d]
        angles[c] = _radians(angle)
        cons[i] = LineInternalAngleConstraint(lines + c, angles + c)
        i += 1
        c += 1

    # Solve
    if solve(pparameters, params_count, cons, cons_count, Rough) != Succsess:
        raise RuntimeError("No valid Solutions were found from this start point")

    """solved_points: List[
        # R joint
        [p1]: (p1_x, p1_y)
        # P or RP joint
        [p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))
    ]
    """
    cdef list solved_points = []
    for i in range(len(vpoints)):
        if vpoints[i].type == VJoint.R:
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
    free(slider_offset)
    free(distances)
    free(angles)
    free(lines)
    free(cons)

    return solved_points
