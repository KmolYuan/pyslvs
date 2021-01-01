# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""PMKS symbolics.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from libc.math cimport M_PI, atan2, hypot
from cpython.object cimport Py_EQ, Py_NE
from numpy import array, zeros, float64 as f64


cdef double distance(double x1, double y1, double x2, double y2) nogil:
    """Distance of two coordinates."""
    return hypot(x2 - x1, y2 - y1)


cdef double slope_angle(double x1, double y1, double x2, double y2) nogil:
    """Slope angle of two coordinates."""
    return atan2(y1 - y2, x1 - x2)


cpdef list get_vlinks(object vpoints):
    """Get VLinks from a list of VPoint `vpoints`."""
    links = {}
    cdef int i
    cdef VPoint vpoint
    for i, vpoint in enumerate(vpoints):
        for name in vpoint.links:
            if name not in links:
                links[name] = {i}
            else:
                links[name].add(i)
    vlinks = []
    for name, points in links.items():
        vlinks.append(VLink(name, "", points))
    return vlinks


@cython.final
cdef class Coord:
    """A data class used to store coordinates."""

    def __cinit__(self, double x, double y):
        self.x = x
        self.y = y

    cpdef double distance(self, Coord p):
        """Return the distance between two coordinates."""
        return distance(self.x, self.y, p.x, p.y)

    cpdef double slope_angle(self, Coord p):
        """Slope angle of two coordinates."""
        return slope_angle(self.x, self.y, p.x, p.y)

    cpdef bint is_nan(self):
        """Return true if the coordinate value is not a number."""
        return self.x != self.x

    def __repr__(self):
        return f"Coord({self.x:.02f}, {self.y:.02f})"


@cython.final
cdef class VPoint:
    """Mechanism expression class."""
    HOLDER = VPoint([], VJoint.R, 0., "", 0., 0.)

    def __cinit__(
        self,
        object links,
        VJoint j_type,
        double angle,
        str color_str,
        double x,
        double y,
        object color_func=None
    ):
        self.set_links(links)
        self.type = j_type
        self.type_str = ('R', 'P', 'RP')[j_type]
        self.angle = angle
        self.color_str = color_str
        if color_func is None:
            self.color = None
        else:
            self.color = color_func(color_str)
        self.x = x
        self.y = y
        self.c = zeros((2, 2), dtype=f64)
        if self.type in {VJoint.P, VJoint.RP}:
            # Slider current coordinates
            # [0]: Current node on slot
            # [1]: Pin
            self.c[0, 0] = self.c[1, 0] = self.x
            self.c[0, 1] = self.c[1, 1] = self.y
        else:
            self.c[0, 0] = self.x
            self.c[0, 1] = self.y
        self.__has_offset = False
        self.__offset = 0

    @staticmethod
    def r_joint(links, x, y):
        """A fast constructor of revolute joints."""
        return VPoint.c_r_joint(links, x, y)

    @staticmethod
    cdef VPoint c_r_joint(object links, double x, double y):
        return VPoint.__new__(VPoint, links, VJoint.R, 0., '', x, y)

    @staticmethod
    def slider_joint(links, type_int, angle, x, y):
        """A fast constructor of slider joints."""
        return VPoint.c_slider_joint(links, type_int, angle, x, y)

    @staticmethod
    cdef VPoint c_slider_joint(object links, VJoint type_int, double angle, double x, double y):
        return VPoint.__new__(VPoint, links, type_int, angle, '', x, y)

    cpdef VPoint copy(self):
        """The copy method of the VPoint object."""
        return self.__copy__()

    @property
    def sx(self):
        """X value of slot coordinate."""
        return self.c[0, 0]

    @property
    def sy(self):
        """Y value of slot coordinate."""
        return self.c[0, 1]

    @property
    def cx(self):
        """X value of current coordinate.
        If it's slider, the pin coordinate will be returned.
        """
        if self.type == VJoint.R:
            return self.c[0, 0]
        else:
            return self.c[1, 0]

    @property
    def cy(self):
        """Y value of current coordinate.
        If it's slider, the pin coordinate will be returned.
        """
        if self.type == VJoint.R:
            return self.c[0, 1]
        else:
            return self.c[1, 1]

    cpdef void set_links(self, object links) except *:
        """The update function of links attribute."""
        self.links = tuple([s for s in links if s])

    cpdef void replace_link(self, str link1, str link2) except *:
        """Replace the value in links attribute."""
        self.set_links([link2 if link == link1 else link for link in self.links])

    cpdef void move(self, object c1, object c2 = None) except *:
        """The update function of current coordinate(s).
        The 2nd placement is the pin coordinate of slider joints.

        If there is only one argument for a slider joint,
        the slot and pin coordinates will be set to the same position.
        """
        cdef double x, y
        x, y = c1
        self.c[0, 0] = x
        self.c[0, 1] = y
        if self.type in {VJoint.P, VJoint.RP}:
            if c2 is not None:
                x, y = c2
            self.c[1, 0] = x
            self.c[1, 1] = y

    cpdef void locate(self, double x, double y) except *:
        """The update function of original coordinate."""
        self.x = x
        self.y = y
        self.c[0, 0] = x
        self.c[0, 1] = y

    cpdef void rotate(self, double angle):
        """The update function of angle attribute."""
        self.angle = angle % 180

    cpdef void set_offset(self, double offset):
        """The update function of slider offset.
        It will also enable offset value after called.
        """
        self.__has_offset = True
        self.__offset = offset

    cpdef void disable_offset(self):
        """Disable offset setting of the joint."""
        self.__has_offset = False

    cpdef bint is_slider(self):
        """Return true for slider type."""
        return self.type in {VJoint.P, VJoint.RP}

    cpdef double distance(self, VPoint p):
        """Return the distance between two VPoint objects."""
        on_links = tuple(set(self.links) & set(p.links))
        cdef double m_x = 0
        cdef double m_y = 0
        cdef double p_x = 0
        cdef double p_y = 0
        if on_links:
            if self.type == VJoint.R or self.links[0] == on_links[0]:
                # self is R joint or at base link
                m_x = self.c[0, 0]
                m_y = self.c[0, 1]
            else:
                # At pin joint
                m_x = self.c[1, 0]
                m_y = self.c[1, 1]
            if p.type == VJoint.R or p.links[0] == on_links[0]:
                # p is R joint or at base link
                p_x = p.c[0, 0]
                p_y = p.c[0, 1]
            else:
                # At pin joint
                p_x = p.c[1, 0]
                p_y = p.c[1, 1]
        else:
            m_x = self.c[0, 0]
            m_y = self.c[0, 1]
            p_x = p.c[0, 0]
            p_y = p.c[0, 1]
        return hypot(p_x - m_x, p_y - m_y)

    cpdef bint has_offset(self):
        """Return true if the offset setting is enabled."""
        return self.__has_offset

    cpdef double offset(self):
        """Return the offset constraint value of the joint."""
        return self.__offset

    cpdef double true_offset(self):
        """Return the current offset value of the joint."""
        return hypot(self.c[1, 0] - self.c[0, 0], self.c[1, 1] - self.c[0, 1])

    cpdef double slope_angle(self, VPoint p, int num1 = 2, int num2 = 2):
        """Return the value `hypot(p_x - m_x, p_y - m_y)`,
        where `m_x`, `m_y` is the value of the joint,
        and `p_x`, `p_y` is the value of `p`.

        The option `num1` and `num2` is the position of current coordinate
        attribute.
        """
        cdef double x1, y1, x2, y2
        if num1 > 1:
            x2, y2 = self.x, self.y
        else:
            x2 = self.c[num2, 0]
            y2 = self.c[num2, 1]
        if num2 > 1:
            x1, y1 = p.x, p.y
        else:
            x1 = p.c[num2, 0]
            y1 = p.c[num2, 1]
        return slope_angle(x1, y1, x2, y2) / M_PI * 180

    cpdef Coord link_pos(self, str link):
        """Return the position for the vlink."""
        cdef size_t ind
        if self.type == VJoint.R or self.is_slot_link(link):
            ind = 0
        else:
            ind = 1
        return self.to_coord(ind)

    cpdef bint grounded(self):
        """Return true if the joint pin is connected to ground link."""
        if self.type == VJoint.R:
            return VLink.FRAME in self.links
        elif self.type in {VJoint.P, VJoint.RP}:
            if self.links:
                return self.is_slot_link(VLink.FRAME)
            else:
                return False

    cpdef bint pin_grounded(self):
        """Return true if the point is at the same link."""
        return VLink.FRAME in self.links[1:]

    cpdef bint same_link(self, VPoint p):
        """Return true if the point is at the same link."""
        return bool(set(self.links) & set(p.links))

    cpdef bint no_link(self):
        """Return true if there is no any link in links attribute."""
        return not self.links

    cpdef bint is_slot_link(self, str link):
        """Return true if the slot is on the link `link_name`."""
        if self.type == VJoint.R:
            return False
        if self.links:
            return link == self.links[0]
        else:
            return False

    cpdef str expr(self):
        """Return the literal mechanism expression of the joint."""
        if self.type != VJoint.R:
            type_text = f"{self.type_str}, A[{self.angle}]"
        else:
            type_text = 'R'
        if self.color_str:
            color = f", color[{self.color_str}]"
        else:
            color = ""
        links_text = ", ".join([name for name in self.links])
        x_text = f"{self.x:.4f}".rstrip('0').rstrip('.')
        y_text = f"{self.y:.4f}".rstrip('0').rstrip('.')
        return f"J[{type_text}{color}, P[{x_text}, {y_text}], L[{links_text}]]"

    cpdef Coord to_coord(self, size_t ind):
        """Obtain coordinate by Coord object."""
        return Coord.__new__(Coord, self.c[ind, 0], self.c[ind, 1])

    def __copy__(self):
        cdef VPoint vpoint = VPoint.__new__(
            VPoint,
            self.links,
            self.type,
            self.angle,
            self.color_str,
            self.x,
            self.y
        )
        vpoint.c[:] = self.c[:]
        return vpoint

    def __richcmp__(self, rhs, int op):
        cdef VPoint other = rhs
        cdef bint different = (
            self.links != other.links or
            self.c[0, 0] != other.c[0, 0] or
            self.c[0, 1] != other.c[0, 1] or
            self.c[1, 1] != other.c[1, 1] or
            self.c[1, 1] != other.c[1, 1] or
            self.type != other.type or
            self.x != other.x or
            self.y != other.y or
            self.angle != other.angle
        )
        if op == Py_EQ:
            return not different
        elif op == Py_NE:
            return different
        else:
            raise TypeError(
                f"'{op}' not support between instances of "
                f"{type(self)} and {type(other)}"
            )

    def __getitem__(self, ind):
        cdef int i = ind
        if self.type == VJoint.R:
            return self.c[0, i]
        else:
            return self.c[1, i]

    def __repr__(self):
        return f"VPoint({self.links}, {int(self.type)}, {self.angle}, {array(self.c).tolist()})"


@cython.final
cdef class VLink:
    """Mechanism expression class in link's view."""
    HOLDER = VLink("", "", [])
    FRAME = 'ground'

    def __cinit__(
        self,
        str name,
        str color_str,
        object points,
        object color_func=None
    ):
        self.name = name
        self.color_str = color_str
        if color_func is None:
            self.color = None
        else:
            self.color = color_func(color_str)
        self.points = array(list(points), dtype=int)

    cpdef void set_points(self, object points) except *:
        """The update function of points attribute."""
        self.points = array(list(points), dtype=int)

    cpdef Coord[:] points_pos(self, object vpoints) except *:
        """Get link positions from a VPoint list."""
        coords = []
        cdef int i
        cdef Coord c
        cdef VPoint vpoint
        for i in self.points:
            vpoint = vpoints[i]
            coords.append(vpoint.link_pos(self.name))
        return array(coords, dtype=object)

    def __contains__(self, point: int):
        return point in self.points

    def __repr__(self):
        return f"VLink('{self.name}', {tuple(self.points)}, color_qt)"
