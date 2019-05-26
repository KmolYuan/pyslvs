# -*- coding: utf-8 -*-
# cython: language_level=3, embedsignature=True, cdivision=True

"""Wrapper of BFGS algorithm.

Note of Pointer:
+ In Cython, pointer is more convenient then array.
    Because we can not "new" them or using "const" decorator on size_t.
+ There is NO pointer's "get value" operator in Cython,
    please use "index" operator.
+ Pointers can be plus with C's Integer, but not Python's.
    So please copy or declare to C's Integer.

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport M_PI, cos, sin
from libcpp.list cimport list as clist
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair
from pyslvs.sketch_solve cimport (
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
from pyslvs.expression cimport (
    get_vlinks,
    VJoint,
    VPoint,
    VLink,
    Coordinate,
)


ctypedef fused T:
    double
    Line


cdef inline T *end_ptr(clist[T] *t_list):
    """Get last pointer."""
    return &t_list.back()


cdef inline void _sort_pairs(dict data_dict):
    """Sort the pairs in data_dict."""
    cdef object k
    for k in data_dict:
        if type(k) is tuple:
            data_dict[frozenset(k)] = data_dict.pop(k)


cdef inline double _radians(double degree):
    """Degrees to radians."""
    return degree / 180 * M_PI


cpdef list vpoint_solving(
    object vpoints_,
    dict inputs = None,
    dict data_dict = None
):
    """Solving function from vpoint list.

    + vpoints: Sequence[VPoint]
    + inputs: {(b0, d0): a0, (b1, d1): a1, ...}

    Known coordinates import from data_dict.
    + data_dict: {0: Coordinate(10.0, 20.0), ..., (0, 2): 30.0, ...}
    """
    # Sequences.
    cdef list vpoints = list(vpoints_)
    cdef dict vlinks = {vlink.name: vlink for vlink in get_vlinks(vpoints)}
    if inputs is None:
        inputs = {}
    if data_dict is None:
        data_dict = {}

    _sort_pairs(data_dict)

    cdef clist[double] params
    cdef clist[double] constants
    cdef vector[Point] points
    cdef cmap[int, int] sliders
    cdef vector[Point] slider_bases
    cdef vector[Point] slider_slots
    cdef clist[Line] slider_lines

    # Point parameters
    cdef int i
    cdef double x, y
    cdef double *tmp_ptr
    cdef VPoint vpoint
    cdef Coordinate coord
    for i, vpoint in enumerate(vpoints):
        if vpoint.no_link():
            x, y = vpoint.c[0]
            constants.push_back(x)
            tmp_ptr = end_ptr(&constants)
            constants.push_back(y)
            points.push_back([tmp_ptr, end_ptr(&constants)])
            continue

        if vpoint.grounded():
            if i in data_dict:
                # Known coordinates.
                coord = data_dict[i]
                constants.push_back(coord.x)
                tmp_ptr = end_ptr(&constants)
                constants.push_back(coord.y)
                points.push_back([tmp_ptr, end_ptr(&constants)])
                continue

            x, y = vpoint.c[0]
            constants.push_back(x)
            tmp_ptr = end_ptr(&constants)
            constants.push_back(y)
            if vpoint.type in {VJoint.P, VJoint.RP}:
                sliders[i] = <int>slider_bases.size()
                # Base point (slot) is fixed.
                slider_bases.push_back([tmp_ptr, end_ptr(&constants)])
                # Slot point (slot) is movable.
                params.push_back(x + cos(vpoint.angle))
                tmp_ptr = end_ptr(&params)
                params.push_back(y + sin(vpoint.angle))
                slider_slots.push_back([tmp_ptr, end_ptr(&params)])
                # Pin is movable.
                x, y = vpoint.c[1]
                if vpoint.has_offset() and vpoint.true_offset() <= 0.1:
                    if vpoint.offset() > 0:
                        x += 0.1
                        y += 0.1
                    else:
                        x -= 0.1
                        y -= 0.1
                params.push_back(x)
                tmp_ptr = end_ptr(&params)
                params.push_back(y)
                points.push_back([tmp_ptr, end_ptr(&params)])
            else:
                points.push_back([tmp_ptr, end_ptr(&constants)])
            continue

        if i in data_dict:
            # Known coordinates.
            coord = data_dict[i]
            constants.push_back(coord.x)
            tmp_ptr = end_ptr(&constants)
            constants.push_back(coord.y)
            points.push_back([tmp_ptr, end_ptr(&constants)])
            continue

        x, y = vpoint.c[0]
        params.push_back(x)
        tmp_ptr = end_ptr(&params)
        params.push_back(y)
        if vpoint.type in {VJoint.P, VJoint.RP}:
            sliders[i] = <int>slider_bases.size()
            # Base point (slot) is movable.
            slider_bases.push_back([tmp_ptr, end_ptr(&params)])
            # Slot point (slot) is movable.
            params.push_back(x + cos(vpoint.angle))
            tmp_ptr = end_ptr(&params)
            params.push_back(y + sin(vpoint.angle))
            slider_slots.push_back([tmp_ptr, end_ptr(&params)])
            if vpoint.pin_grounded():
                # Pin is fixed.
                x, y = vpoint.c[1]
                constants.push_back(x)
                tmp_ptr = end_ptr(&constants)
                constants.push_back(y)
                points.push_back([tmp_ptr, end_ptr(&constants)])
            else:
                # Pin is movable.
                x, y = vpoint.c[1]
                if vpoint.has_offset() and vpoint.true_offset() <= 0.1:
                    if vpoint.offset() > 0:
                        x += 0.1
                        y += 0.1
                    else:
                        x -= 0.1
                        y -= 0.1
                params.push_back(x)
                tmp_ptr = end_ptr(&params)
                params.push_back(y)
                points.push_back([tmp_ptr, end_ptr(&params)])
            continue

        # Point is movable.
        points.push_back([tmp_ptr, end_ptr(&params)])

    cdef clist[Constraint] cons_list

    # Link constraints
    cdef int a, b, c, d
    cdef frozenset frozen_pair
    cdef VPoint vp1, vp2
    cdef Point *p1
    cdef Point *p2
    cdef VLink vlink
    for vlink in vlinks.values():
        if len(vlink.points) < 2:
            continue
        if vlink.name == 'ground':
            continue

        a = vlink.points[0]
        b = vlink.points[1]
        if (a not in data_dict) or (b not in data_dict):
            vp1 = vpoints[a]
            vp2 = vpoints[b]
            if a not in data_dict and vp1.is_slot_link(vlink.name):
                p1 = &slider_bases[sliders[a]]
            else:
                p1 = &points[a]

            if b not in data_dict and vp2.is_slot_link(vlink.name):
                p2 = &slider_bases[sliders[b]]
            else:
                p2 = &points[b]

            frozen_pair = frozenset({a, b})
            if frozen_pair in data_dict:
                constants.push_back(data_dict[frozen_pair])
            else:
                constants.push_back(vp1.distance(vp2))

            cons_list.push_back(P2PDistanceConstraint(p1, p2, end_ptr(&constants)))

        for c in vlink.points[2:]:
            if c in data_dict:
                # Known coordinate.
                continue
            for d in (a, b):
                vp1 = vpoints[c]
                vp2 = vpoints[d]
                if vp1.is_slot_link(vlink.name):
                    p1 = &slider_bases[sliders[c]]
                else:
                    p1 = &points[c]

                if d not in data_dict and vp2.is_slot_link(vlink.name):
                    p2 = &slider_bases[sliders[d]]
                else:
                    p2 = &points[d]

                frozen_pair = frozenset({c, d})
                if frozen_pair in data_dict:
                    constants.push_back(data_dict[frozen_pair])
                else:
                    constants.push_back(vp1.distance(vp2))

                cons_list.push_back(P2PDistanceConstraint(p1, p2, end_ptr(&constants)))

    # Slider constraints
    cdef str name
    cdef Line *slider_slot
    cdef pair[int, int] slider
    for slider in sliders:
        a = slider.first
        b = slider.second
        # Base point
        vp1 = vpoints[a]
        p1 = &points[a]
        # Base slot
        slider_lines.push_back([&slider_bases[b], &slider_slots[b]])
        slider_slot = end_ptr(&slider_lines)
        if vp1.grounded():
            # Slot is grounded.
            constants.push_back(_radians(vp1.angle))
            cons_list.push_back(LineInternalAngleConstraint(slider_slot, end_ptr(&constants)))
            cons_list.push_back(PointOnLineConstraint(p1, slider_slot))
            if vp1.has_offset():
                p2 = &slider_bases[b]
                if vp1.offset():
                    constants.push_back(vp1.offset())
                    cons_list.push_back(P2PDistanceConstraint(p2, p1, end_ptr(&constants)))
                else:
                    cons_list.push_back(PointOnPointConstraint(p2, p1))
        else:
            # Slider between links.
            for name in vp1.links[:1]:
                vlink = vlinks[name]
                # A base link friend.
                c = vlink.points[0]
                if c == a:
                    if len(vlink.points) < 2:
                        # If no any friend.
                        continue
                    c = vlink.points[1]

                vp2 = vpoints[c]
                if vp2.is_slot_link(vlink.name):
                    # c is a slider, and it is be connected with slot link.
                    p2 = &slider_bases[sliders[c]]
                else:
                    # c is a R joint or it is not connected with slot link.
                    p2 = &points[c]
                slider_lines.push_back([&slider_bases[b], p2])
                constants.push_back(_radians(vp1.slope_angle(vp2) - vp1.angle))
                cons_list.push_back(InternalAngleConstraint(
                    slider_slot,
                    end_ptr(&slider_lines),
                    end_ptr(&constants)
                ))
                cons_list.push_back(PointOnLineConstraint(p1, slider_slot))

                if vp1.has_offset():
                    p2 = &slider_bases[b]
                    if vp1.offset():
                        constants.push_back(vp1.offset())
                        cons_list.push_back(P2PDistanceConstraint(p2, p1, end_ptr(&constants)))
                    else:
                        cons_list.push_back(PointOnPointConstraint(p2, p1))

        if vp1.type != VJoint.P:
            continue

        for name in vp1.links[1:]:
            vlink = vlinks[name]
            # A base link friend.
            c = vlink.points[0]
            if c == a:
                if len(vlink.points) < 2:
                    # If no any friend.
                    continue
                c = vlink.points[1]

            vp2 = vpoints[c]
            if vp2.is_slot_link(vlink.name):
                # c is a slider, and it is be connected with slot link.
                p2 = &slider_bases[sliders[c]]
            else:
                # c is a R joint or it is not connected with slot link.
                p2 = &points[c]
            slider_lines.push_back([p1, p2])
            constants.push_back(_radians(vp1.slope_angle(vp2) - vp1.angle))
            cons_list.push_back(InternalAngleConstraint(
                slider_slot,
                end_ptr(&slider_lines),
                end_ptr(&constants)
            ))

    # Angle constraints
    cdef clist[Line] handles
    cdef double angle
    for (b, d), angle in inputs.items():
        if b == d:
            continue
        handles.push_back([&points[b], &points[d]])
        constants.push_back(_radians(angle))
        cons_list.push_back(LineInternalAngleConstraint(
            end_ptr(&handles),
            end_ptr(&constants)
        ))

    # Pointer of parameters
    cdef int params_count = <int>params.size()
    cdef double **params_ptr = <double **>PyMem_Malloc(sizeof(double *) * params_count)
    cdef clist[double].iterator it = params.begin()
    for i in range(params_count):
        params_ptr[i] = &cython.operator.dereference(cython.operator.postincrement(it))

    # Pointer of constraints
    cdef int cons_count = <int>cons_list.size()
    cdef Constraint *cons = <Constraint *>PyMem_Malloc(sizeof(Constraint) * cons_count)
    i = 0
    cdef Constraint con
    for con in cons_list:
        cons[i] = con
        i += 1

    # Solve
    cdef int flag = solve(params_ptr, params_count, cons, cons_count, Rough)

    cdef list solved_points
    if flag == Succsess:
        solved_points = []
        for i, vpoint in enumerate(vpoints):
            if vpoint.type == VJoint.R:
                solved_points.append((points[i].x[0], points[i].y[0]))
            else:
                solved_points.append((
                    (slider_bases[sliders[i]].x[0], slider_bases[sliders[i]].y[0]),
                    (points[i].x[0], points[i].y[0])
                ))

    PyMem_Free(params_ptr)
    PyMem_Free(cons)

    if flag == Succsess:
        return solved_points
    else:
        raise ValueError("no valid solutions were found from initialed values")
