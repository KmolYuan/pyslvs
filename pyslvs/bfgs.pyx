# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Wrapper of BFGS algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from typing import Sequence, Set
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport M_PI, cos, sin
from libcpp.pair cimport pair
from .sketch_solve cimport (
    point_on_point,
    p2p_distance,
    internal_angle,
    point_on_line,
    line_internal_angle,
    solve,
)
from .expression cimport get_vlinks, VJoint, VPoint, Coord


cdef inline double *de_refer_post_inc(clist[double].iterator &it):
    """Implement &(*it++) in C++."""
    return &cython.operator.dereference(cython.operator.postincrement(it))


cdef inline void _sort_pairs(dict data_dict):
    """Sort the pairs in data_dict."""
    for k in data_dict:
        if isinstance(k, (Sequence, Set)) and len(k) == 2:
            data_dict[frozenset(k)] = data_dict.pop(k)


cdef inline double _radians(double degree):
    """Degrees to radians."""
    return degree / 180 * M_PI


cdef class SolverSystem:
    """Sketch Solve solver.

    !!! note
        The object attributes of such type are unable to access.
    """

    def __cinit__(self, object vpoints_, dict inputs = None, dict data_dict = None):
        self.vpoints = list(vpoints_)
        self.vlinks = {vlink.name: vlink for vlink in get_vlinks(self.vpoints)}
        self.inputs = {} if inputs is None else inputs
        self.data_dict = {} if data_dict is None else data_dict
        _sort_pairs(self.data_dict)
        self.build_expression()

    cpdef bint same_points(self, object vpoints_):
        """Return true if two expressions are same."""
        cdef int i
        cdef VPoint p1, p2
        for i, p1 in enumerate(vpoints_):
            p2 = self.vpoints[i]
            if p1.links != p2.links:
                return False
        return True

    cpdef frozenset show_inputs(self):
        """Show the current input pairs keys from original constructor."""
        return frozenset(self.inputs)

    cpdef frozenset show_data(self):
        """Show the current keys of `data_dict` parameter from original
        constructor.
        """
        return frozenset(self.data_dict)

    cdef Point *point_ptr(self, int i, VLink vlink):
        """Pick Point pointers."""
        cdef VPoint vp = self.vpoints[i]
        if i not in self.data_dict and vp.is_slot_link(vlink.name):
            return &self.slider_bases[self.sliders[i]]
        else:
            return &self.points[i]

    cdef void build_expression(self):
        """Build the expression for solver at first time."""
        # Point parameters
        cdef int i
        cdef double x, y
        cdef double *tmp_ptr
        cdef VPoint vpoint
        cdef Coord coord
        for i, vpoint in enumerate(self.vpoints):
            if vpoint.no_link():
                x = vpoint.c[0, 0]
                y = vpoint.c[0, 1]
                self.constants.push_back(x)
                tmp_ptr = &self.constants.back()
                self.constants.push_back(y)
                self.points.push_back([tmp_ptr, &self.constants.back()])
                continue
            if vpoint.grounded():
                if self.check_known(i):
                    continue
                x = vpoint.c[0, 0]
                y = vpoint.c[0, 1]
                self.constants.push_back(x)
                tmp_ptr = &self.constants.back()
                self.constants.push_back(y)
                if vpoint.type in {VJoint.P, VJoint.RP}:
                    self.sliders[i] = <int>self.slider_bases.size()
                    # Base point (slot) is fixed
                    self.slider_bases.push_back([tmp_ptr, &self.constants.back()])
                    # Slot point (slot) is movable
                    self.params.push_back(x + cos(vpoint.angle))
                    tmp_ptr = &self.params.back()
                    self.params.push_back(y + sin(vpoint.angle))
                    self.slider_slots.push_back([tmp_ptr, &self.params.back()])
                    # Pin is movable
                    x = vpoint.c[1, 0]
                    y = vpoint.c[1, 1]
                    if vpoint.has_offset() and vpoint.true_offset() <= 0.1:
                        if vpoint.offset() > 0:
                            x += 0.1
                            y += 0.1
                        else:
                            x -= 0.1
                            y -= 0.1
                    self.params.push_back(x)
                    tmp_ptr = &self.params.back()
                    self.params.push_back(y)
                    self.points.push_back([tmp_ptr, &self.params.back()])
                else:
                    self.points.push_back([tmp_ptr, &self.constants.back()])
                continue
            if self.check_known(i):
                continue
            x = vpoint.c[0, 0]
            y = vpoint.c[0, 1]
            self.params.push_back(x)
            tmp_ptr = &self.params.back()
            self.params.push_back(y)
            if vpoint.type in {VJoint.P, VJoint.RP}:
                self.sliders[i] = <int>self.slider_bases.size()
                # Base point (slot) is movable
                self.slider_bases.push_back([tmp_ptr, &self.params.back()])
                # Slot point (slot) is movable
                self.params.push_back(x + cos(vpoint.angle))
                tmp_ptr = &self.params.back()
                self.params.push_back(y + sin(vpoint.angle))
                self.slider_slots.push_back([tmp_ptr, &self.params.back()])
                if vpoint.pin_grounded():
                    # Pin is fixed
                    x = vpoint.c[1, 0]
                    y = vpoint.c[1, 1]
                    self.constants.push_back(x)
                    tmp_ptr = &self.constants.back()
                    self.constants.push_back(y)
                    self.points.push_back([tmp_ptr, &self.constants.back()])
                else:
                    # Pin is movable
                    x = vpoint.c[1, 0]
                    y = vpoint.c[1, 1]
                    if vpoint.has_offset() and vpoint.true_offset() <= 0.1:
                        if vpoint.offset() > 0:
                            x += 0.1
                            y += 0.1
                        else:
                            x -= 0.1
                            y -= 0.1
                    self.params.push_back(x)
                    tmp_ptr = &self.params.back()
                    self.params.push_back(y)
                    self.points.push_back([tmp_ptr, &self.params.back()])
                continue
            # Point is movable
            self.points.push_back([tmp_ptr, &self.params.back()])
        # Link constraints
        # (automatic fill up the link length options of data keys)
        cdef int a, b, c, d
        cdef VPoint vp1, vp2
        cdef Point *p1
        cdef Point *p2
        cdef VLink vlink
        for vlink in self.vlinks.values():
            if len(vlink.points) < 2 or vlink.name == VLink.FRAME:
                continue
            a = vlink.points[0]
            b = vlink.points[1]
            if (a not in self.data_dict) or (b not in self.data_dict):
                vp1 = self.vpoints[a]
                vp2 = self.vpoints[b]
                p1 = self.point_ptr(a, vlink)
                p2 = self.point_ptr(b, vlink)
                frozen_pair = frozenset({a, b})
                if frozen_pair in self.data_dict:
                    x = self.data_dict[frozen_pair]
                else:
                    x = vp1.distance(vp2)
                    self.data_dict[frozen_pair] = x
                self.data_values.push_back(x)
                self.cons_list.push_back(p2p_distance(p1, p2, &self.data_values.back()))
            for c in vlink.points[2:]:
                if c in self.data_dict:
                    # Known coordinate
                    continue
                for d in (a, b):
                    vp1 = self.vpoints[c]
                    vp2 = self.vpoints[d]
                    p1 = self.point_ptr(c, vlink)
                    p2 = self.point_ptr(d, vlink)
                    frozen_pair = frozenset({c, d})
                    if frozen_pair in self.data_dict:
                        x = self.data_dict[frozen_pair]
                    else:
                        x = vp1.distance(vp2)
                        self.data_dict[frozen_pair] = x
                    self.data_values.push_back(x)
                    self.cons_list.push_back(p2p_distance(p1, p2, &self.data_values.back()))
        # Slider constraints
        cdef Line *slider_slot
        cdef pair[int, int] slider
        for slider in self.sliders:
            a = slider.first
            b = slider.second
            # Base point
            vp1 = self.vpoints[a]
            p1 = &self.points[a]
            # Base slot
            self.slider_lines.push_back([&self.slider_bases[b], &self.slider_slots[b]])
            slider_slot = &self.slider_lines.back()
            if vp1.grounded():
                # Slot is grounded
                self.constants.push_back(_radians(vp1.angle))
                self.cons_list.push_back(line_internal_angle(slider_slot, &self.constants.back()))
                self.cons_list.push_back(point_on_line(p1, slider_slot))
                if vp1.has_offset():
                    p2 = &self.slider_bases[b]
                    if vp1.offset():
                        self.constants.push_back(vp1.offset())
                        self.cons_list.push_back(p2p_distance(p2, p1, &self.constants.back()))
                    else:
                        self.cons_list.push_back(point_on_point(p2, p1))
            else:
                # Slider between links
                for name in vp1.links[:1]:
                    vlink = self.vlinks[name]
                    # A base link friend
                    c = vlink.points[0]
                    if c == a:
                        if len(vlink.points) < 2:
                            # If no any friend
                            continue
                        c = vlink.points[1]
                    vp2 = self.vpoints[c]
                    if vp2.is_slot_link(vlink.name):
                        # c is a slider, and it is be connected with slot link
                        p2 = &self.slider_bases[self.sliders[c]]
                    else:
                        # c is a R joint or it is not connected with slot link
                        p2 = &self.points[c]
                    self.slider_lines.push_back([&self.slider_bases[b], p2])
                    self.constants.push_back(_radians(vp1.slope_angle(vp2) - vp1.angle))
                    self.cons_list.push_back(internal_angle(
                        slider_slot,
                        &self.slider_lines.back(),
                        &self.constants.back()
                    ))
                    self.cons_list.push_back(point_on_line(p1, slider_slot))
                    if vp1.has_offset():
                        p2 = &self.slider_bases[b]
                        if vp1.offset():
                            self.constants.push_back(vp1.offset())
                            self.cons_list.push_back(p2p_distance(p2, p1, &self.constants.back()))
                        else:
                            self.cons_list.push_back(point_on_point(p2, p1))
            if vp1.type != VJoint.P:
                continue
            for name in vp1.links[1:]:
                vlink = self.vlinks[name]
                # A base link friend
                c = vlink.points[0]
                if c == a:
                    if len(vlink.points) < 2:
                        # If no any friend
                        continue
                    c = vlink.points[1]
                vp2 = self.vpoints[c]
                if vp2.is_slot_link(vlink.name):
                    # c is a slider, and it is be connected with slot link
                    p2 = &self.slider_bases[self.sliders[c]]
                else:
                    # c is a R joint or it is not connected with slot link
                    p2 = &self.points[c]
                self.slider_lines.push_back([p1, p2])
                self.constants.push_back(_radians(vp1.slope_angle(vp2) - vp1.angle))
                self.cons_list.push_back(internal_angle(
                    slider_slot,
                    &self.slider_lines.back(),
                    &self.constants.back()
                ))
        # Angle constraints
        cdef double angle
        for (b, d), angle in self.inputs.items():
            if b == d:
                continue
            self.handles.push_back([&self.points[b], &self.points[d]])
            self.inputs_angle.push_back(_radians(angle))
            self.cons_list.push_back(line_internal_angle(
                &self.handles.back(),
                &self.inputs_angle.back()
            ))

    cdef bint check_known(self, int i):
        """Check known coordinates."""
        if i not in self.data_dict:
            return False
        coord = self.data_dict[i]
        self.data_values.push_back(coord.x)
        tmp_ptr = &self.data_values.back()
        self.data_values.push_back(coord.y)
        self.points.push_back([tmp_ptr, &self.data_values.back()])
        return True

    cpdef void set_inputs(self, dict inputs):
        """Set the values of `inputs` parameter from original constructor.
        Two groups of `dict` keys must be the same or subset.
        """
        if self.inputs is None or inputs is None:
            raise ValueError(f"do not accept modifications")
        if not self.show_inputs() >= set(inputs):
            raise ValueError(f"format must be {set(self.inputs)}, not {set(inputs)}")
        self.inputs.update(inputs)
        # Set values
        cdef int b, d
        cdef double angle
        cdef double *handle
        cdef clist[double].iterator it = self.inputs_angle.begin()
        for (b, d), angle in self.inputs.items():
            if b == d:
                continue
            handle = de_refer_post_inc(it)
            handle[0] = _radians(angle)

    cpdef void set_data(self, dict data_dict):
        """Set the values of `data_dict` parameter from original constructor.
        Two groups of `dict` keys must be the same or subset.
        """
        if self.data_dict is None or data_dict is None:
            raise ValueError(f"do not accept modifications")
        _sort_pairs(data_dict)
        if not self.show_data() >= set(data_dict):
            raise ValueError(f"format must be {set(self.data_dict)}, not {set(data_dict)}")
        self.data_dict.update(data_dict)
        cdef size_t n = 0
        # Set values
        cdef int i
        cdef double *handle
        cdef VPoint vpoint
        cdef Coord coord
        cdef clist[double].iterator it = self.data_values.begin()
        for i, vpoint in enumerate(self.vpoints):
            if vpoint.grounded():
                if i in self.data_dict:
                    # Known coordinates
                    coord = self.data_dict[i]
                    handle = de_refer_post_inc(it)
                    handle[0] = coord.x
                    handle = de_refer_post_inc(it)
                    handle[0] = coord.y
            if i in self.data_dict:
                # Known coordinates
                coord = self.data_dict[i]
                handle = de_refer_post_inc(it)
                handle[0] = coord.x
                handle = de_refer_post_inc(it)
                handle[0] = coord.y
        cdef int a, b, c, d
        cdef VLink vlink
        for vlink in self.vlinks.values():
            if len(vlink.points) < 2:
                continue
            if vlink.name == VLink.FRAME:
                continue
            a = vlink.points[0]
            b = vlink.points[1]
            if (a not in self.data_dict) or (b not in self.data_dict):
                handle = de_refer_post_inc(it)
                handle[0] = self.data_dict[frozenset({a, b})]
            for c in vlink.points[2:]:
                if c in self.data_dict:
                    # Known coordinate
                    continue
                for d in (a, b):
                    handle = de_refer_post_inc(it)
                    handle[0] = self.data_dict[frozenset({c, d})]

    cpdef list solve(self):
        """Solve the conditions and return the result, raise ValueError if
        not succeeded.
        The joint position will returned by its index correspondingly.

        + Revolute joints: Tuple[float, float]
        + Slider joints: Tuple[Tuple[float, float], Tuple[float, float]]
        """
        # Pointer of parameters
        cdef size_t params_count = <int>self.params.size()
        cdef double **params_ptr = <double **>PyMem_Malloc(sizeof(double *) * params_count)
        cdef clist[double].iterator it = self.params.begin()
        cdef size_t i
        for i in range(params_count):
            params_ptr[i] = de_refer_post_inc(it)
        # Pointer of constraints
        cdef size_t cons_count = <int>self.cons_list.size()
        cdef Constraint *cons = <Constraint *>PyMem_Malloc(sizeof(Constraint) * cons_count)
        i = 0
        cdef Constraint con
        for con in self.cons_list:
            cons[i] = con
            i += 1
        # Solve
        cdef bint flag = solve(params_ptr, params_count, cons, cons_count, False)
        cdef VPoint vp
        cdef Coord c
        cdef Point p1, p2
        if flag:
            solved_points = []
            for i, vp in enumerate(self.vpoints):
                if i in self.data_dict or vp.no_link():
                    if vp.no_link():
                        c = Coord.__new__(Coord, vp.c[0, 0], vp.c[0, 1])
                    else:
                        c = self.data_dict[i]
                    if vp.type == VJoint.R:
                        solved_points.append((c.x, c.y))
                    else:
                        solved_points.append(((c.x, c.y), (c.x, c.y)))
                    continue
                p1 = self.points[i]
                if vp.type == VJoint.R:
                    solved_points.append((p1.x[0], p1.y[0]))
                else:
                    p2 = self.slider_bases[self.sliders[i]]
                    solved_points.append(((p2.x[0], p2.y[0]),
                                          (p1.x[0], p1.y[0])))
        PyMem_Free(params_ptr)
        PyMem_Free(cons)
        if flag:
            return solved_points
        else:
            raise ValueError("no valid solutions were found from initialed values")
