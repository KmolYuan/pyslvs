# -*- coding: utf-8 -*-
# cython: language_level=3

"""PMKS simbolics."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from libc.math cimport (
    M_PI,
    atan2,
    hypot,
)
from numpy cimport ndarray
from numpy import object as np_object


cdef class VPoint:
    
    """Symbol of joints."""
    
    def __cinit__(self,
        links: str,
        type_int: int,
        angle: double,
        color_str: str,
        x: double,
        y: double,
        color_func: object = None
    ):
        cdef list tmp_list = []
        cdef str name
        links = links.replace(" ", '')
        for name in links.split(','):
            if not name:
                continue
            tmp_list.append(name)
        self.links = tuple(tmp_list)
        self.type = type_int
        self.typeSTR = ('R', 'P', 'RP')[type_int]
        self.angle = angle
        self.colorSTR = color_str
        if color_func:
            self.color = color_func(color_str)
        self.x = x
        self.y = y
        self.c = ndarray(2, dtype=np_object)
        if (self.type == 1) or (self.type == 2):
            """Slider current coordinates.
            
            + [0]: Current node on slot.
            + [1]: Pin.
            """
            self.c[0] = (self.x, self.y)
            self.c[1] = (self.x, self.y)
        else:
            self.c[0] = (self.x, self.y)
    
    @property
    def cx(self):
        """X value of frist current coordinate."""
        return self.c[0][0]
    
    @property
    def cy(self):
        """Y value of frist current coordinate."""
        return self.c[0][1]
    
    cpdef void move(self, tuple c1, tuple c2 = None):
        """Change coordinates of this point."""
        self.c[0] = c1
        if self.type != 0:
            self.c[1] = c2 if c2 else c1
    
    cpdef void rotate(self, double angle):
        """Change the angle of slider slot by degrees."""
        self.angle = angle % 180
    
    cpdef double distance(self, VPoint p):
        """Distance between two VPoint."""
        return hypot(p.x - self.x, p.y - self.y)
    
    cpdef double slope_angle(self, VPoint p, int num1 = 2, int num2 = 2):
        """Angle between horizontal line and two point.
        
        num1: me.
        num2: other side.
        [0]: base (slot) link.
        [1]: pin link.
        """
        cdef double x1, y1, x2, y2
        if num1 > 1:
            x2, y2 = self.x, self.y
        else:
            x2, y2 = self.c[num2]
        if num2 > 1:
            x1, y1 = p.x, p.y
        else:
            x1, y1 = p.c[num2]
        return atan2(y1 - y2, x1 - x2) / M_PI * 180
    
    cpdef bool grounded(self):
        """Return True if the joint is connect with the ground."""
        if self.type == 0:
            return 'ground' in self.links
        elif self.type in {1, 2}:
            if self.links:
                return 'ground' == self.links[0]
            else:
                return False
    
    cpdef bool is_slot_link(self, str link_name):
        """Return True if the link name is first link."""
        if self.type == 0:
            return False
        if self.links:
            return link_name == self.links[0]
        else:
            return False
    
    @property
    def expr(self):
        """Expression."""
        return "J[{}, color[{}], P[{}], L[{}]]".format(
            f"{self.typeSTR}, A[{self.angle}]" if self.typeSTR != 'R' else 'R',
            self.colorSTR,
            f"{self.x}, {self.y}",
            ", ".join(l for l in self.links)
        )
    
    def __getitem__(self, i: int) -> float:
        """Get coordinate like this:
        
        x, y = VPoint(10, 20)
        """
        if self.type == 0:
            return self.c[0][i]
        else:
            return self.c[1][i]
    
    def __repr__(self):
        """Use to generate script."""
        return f"VPoint({self.links}, {self.type}, {self.angle}, {self.c})"


cdef class VLink:
    
    """Symbol of links."""
    
    cdef readonly str name, colorSTR
    cdef readonly object color
    cdef readonly tuple points
    
    def __cinit__(self,
        str name,
        str color_str,
        tuple points,
        object color_func = None
    ):
        self.name = name
        self.colorSTR = color_str
        if color_func:
            self.color = color_func(color_str)
        self.points = points
    
    def __contains__(self, point: int):
        """Check if point number is in the link."""
        return point in self.points
    
    def __repr__(self):
        """Use to generate script."""
        return f"VLink('{self.name}', {self.points}, colorQt)"
