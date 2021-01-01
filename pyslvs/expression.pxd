# -*- coding: utf-8 -*-
# cython: language_level=3

"""PMKS symbolics.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cdef double distance(double x1, double y1, double x2, double y2) nogil
cdef double slope_angle(double x1, double y1, double x2, double y2) nogil
cpdef list get_vlinks(object vpoints)

cdef struct CCoord:
    double x, y

cdef class Coord:
    # Coordinate data type.
    # Needs to replaced at backend with c-struct type CCoord.
    cdef public double x, y

    cpdef double distance(self, Coord p)
    cpdef double slope_angle(self, Coord p)
    cpdef bint is_nan(self)

cpdef enum VJoint:
    # Joint types.
    # Actually "class VJoint(IntEnum)" in Python but "enum" in C++.
    R  # Rotate pair
    P  # Prismatic pair
    RP  # Rotate and prismatic pair

cdef class VPoint:
    # VPoint(links, type_int, angle, color_str, x, y, color_func=None)

    cdef readonly tuple links
    cdef readonly double[:, :] c
    cdef readonly VJoint type
    cdef readonly tuple color
    cdef readonly str color_str
    cdef readonly str type_str
    cdef readonly double x, y, angle
    cdef double __offset
    cdef bint __has_offset

    @staticmethod
    cdef VPoint c_r_joint(object links, double x, double y)
    @staticmethod
    cdef VPoint c_slider_joint(object links, VJoint type_int, double angle, double x, double y)

    # Copy method
    cpdef VPoint copy(self)

    # Set values
    cpdef void set_links(self, object links) except *
    cpdef void replace_link(self, str link1, str link2) except *
    cpdef void move(self, object c1, object c2 = *) except *
    cpdef void locate(self, double x, double y) except *
    cpdef void rotate(self, double angle)
    cpdef void set_offset(self, double offset)
    cpdef void disable_offset(self)

    # Get or calculate values
    cpdef bint is_slider(self)
    cpdef double distance(self, VPoint p)
    cpdef bint has_offset(self)
    cpdef double offset(self)
    cpdef double true_offset(self)
    cpdef double slope_angle(self, VPoint p, int num1 = *, int num2 = *)
    cpdef Coord link_pos(self, str link)

    # Link operators
    cpdef bint grounded(self)
    cpdef bint pin_grounded(self)
    cpdef bint same_link(self, VPoint p)
    cpdef bint no_link(self)
    cpdef bint is_slot_link(self, str link)

    # Expression
    cpdef str expr(self)
    cpdef Coord to_coord(self, size_t ind)

cdef class VLink:
    # VLink(name, color_str, points, color_func=None)

    cdef readonly str name, color_str
    cdef readonly tuple color
    cdef readonly long[:] points

    # Set values
    cpdef void set_points(self, object points) except *
    cpdef Coord[:] points_pos(self, object vpoints) except *
