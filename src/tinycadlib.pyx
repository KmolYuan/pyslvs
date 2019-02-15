# -*- coding: utf-8 -*-
# cython: language_level=3

"""Tiny CAD library of PMKS simbolic and position analysis.

author: Yuan Chang
copyright: Copyright (C) 2016-2019
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from libc.math cimport (
    M_PI,
    sqrt,
    sin,
    cos,
    atan2,
    hypot,
    isnan,
    NAN,
)
from expression cimport VJoint, VPoint
from bfgs cimport vpoint_solving


@cython.cdivision
cdef inline double radians(double degree) nogil:
    """Deg to rad."""
    return degree / 180 * M_PI


cdef inline double distance(double x1, double y1, double x2, double y2) nogil:
    """Distance of two cartesian coordinates."""
    return hypot(x2 - x1, y2 - y1)


@cython.final
cdef class Coordinate:

    """A class to store the coordinate."""

    def __cinit__(self, x: double, y: double):
        self.x = x
        self.y = y

    cpdef double distance(self, Coordinate p):
        """Distance."""
        return distance(self.x, self.y, p.x, p.y)

    cpdef bint is_nan(self):
        """Test this coordinate is a error-occurred answer."""
        return isnan(self.x)

    def __repr__(self):
        """Debug printing."""
        return f"Coordinate({self.x:.02f}, {self.y:.02f})"


cpdef tuple plap(
    Coordinate c1,
    double d0,
    double a0,
    Coordinate c2 = None,
    bint inverse = False
):
    """Point on circle by angle."""
    cdef double a1 = atan2(c2.y - c1.y, c2.x - c1.x) if c2 is not None else 0
    if inverse:
        return (c1.x + d0 * cos(a1 - a0)), (c1.y + d0 * sin(a1 - a0))
    else:
        return (c1.x + d0 * cos(a1 + a0)), (c1.y + d0 * sin(a1 + a0))


@cython.cdivision
cpdef tuple pllp(
    Coordinate c1,
    double d0,
    double d1,
    Coordinate c2,
    bint inverse = False
):
    """Two intersection points of two circles."""
    cdef double dx = c2.x - c1.x
    cdef double dy = c2.y - c1.y
    cdef double d = c1.distance(c2)

    # No solutions, the circles are separate.
    if d > d0 + d1:
        return NAN, NAN

    # No solutions because one circle is contained within the other.
    if d < abs(d0 - d1):
        return NAN, NAN

    # Circles are coincident and there are an infinite number of solutions.
    if d == 0 and d0 == d1:
        return NAN, NAN
    cdef double a = (d0 * d0 - d1 * d1 + d * d) / (2 * d)
    cdef double h = sqrt(d0 * d0 - a * a)
    cdef double xm = c1.x + a * dx / d
    cdef double ym = c1.y + a * dy / d

    if inverse:
        return (xm + h * dy / d), (ym - h * dx / d)
    else:
        return (xm - h * dy / d), (ym + h * dx / d)


@cython.cdivision
cpdef tuple plpp(
    Coordinate c1,
    double d0,
    Coordinate c2,
    Coordinate c3,
    bint inverse = False
):
    """Two intersection points of a line and a circle."""
    cdef double line_mag = c2.distance(c3)
    cdef double dx = c3.x - c2.x
    cdef double dy = c3.y - c2.y
    cdef double u = ((c1.x - c2.x) * dx + (c1.y - c2.y) * dy) / (line_mag * line_mag)
    cdef Coordinate inter = Coordinate(c2.x + u * dx, c2.y + u * dy)

    # Test distance between point A and intersection.
    cdef double d = c1.distance(inter)
    if d > d0:
        # No intersection.
        return NAN, NAN
    elif d == d0:
        # One intersection point.
        return inter.x, inter.y

    # Two intersection points.
    d = sqrt(d0 * d0 - d * d) / line_mag
    if inverse:
        return (inter.x - dx * d), (inter.y - dy * d)
    else:
        return (inter.x + dx * d), (inter.y + dy * d)


cpdef tuple pxy(Coordinate c1, double x, double y):
    """Using relative cartesian coordinate to get solution."""
    return (c1.x + x), (c1.y + y)


cdef inline str str_between(str s, str front, str back):
    """Get the string that is inside of parenthesis."""
    return s[s.find(front) + 1:s.find(back)]


cdef inline str str_before(str s, str front):
    """Get the string that is front of parenthesis."""
    return s[:s.find(front)]


cpdef void expr_parser(object exprs, dict data_dict):
    """Update data.
    
    + exprs: "PLAP[P0,L0,a0,P1](P2);PLLP[P2,L1,L2,P1](P3);..."
    + data_dict: {'a0':0., 'L1':10., 'A':(30., 40.), ...}
    """
    cdef object expr
    cdef str func, p
    cdef list params
    cdef int params_count
    cdef double x, y
    for expr in exprs:
        # If the mechanism has no any solution.
        if not expr:
            return
        func = expr[0]
        params = []
        for p in expr[1:-1]:
            if p == 'T':
                params.append(True)
            elif p == 'F':
                params.append(False)
            elif type(data_dict[p]) == tuple:
                x, y = data_dict[p]
                params.append(Coordinate(x, y))
            else:
                params.append(data_dict[p])
        params_count = len(params)

        # We should unpack as C++'s way.
        x = NAN
        y = NAN
        if func == 'PLAP':
            if params_count == 3:
                x, y = plap(params[0], params[1], params[2])
            elif params_count == 4:
                x, y = plap(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = plap(params[0], params[1], params[2], params[3], params[4])
        elif func == 'PLLP':
            if params_count == 4:
                x, y = pllp(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = pllp(params[0], params[1], params[2], params[3], params[4])
        elif func == 'PLPP':
            if params_count == 4:
                x, y = plpp(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = plpp(params[0], params[1], params[2], params[3], params[4])
        elif func == 'PXY':
            if params_count == 3:
                x, y = pxy(params[0], params[1], params[2])
        data_dict[expr[-1]] = (x, y)


cpdef str expr_join(object exprs):
    """Use to append a list of symbols into a string."""
    return ';'.join([
        f"{expr[0]}[{','.join(expr[1:-1])}, {expr[-1]}]({expr[-1]})" for expr in exprs
    ])


cpdef tuple expr_parse(str exprs):
    """Parse expression as tuple."""
    exprs = exprs.replace(" ", '')
    cdef list tmp_list = []
    cdef list params = []
    cdef str expr, p
    for expr in exprs.split(';'):
        if not expr:
            return ()
        params.clear()
        params.append(str_before(expr, '['))
        params.extend(str_between(expr, '[', ']').split(','))
        params.append(str_between(expr, '(', ')'))
        tmp_list.append(tuple(params))
    return tuple(tmp_list)


cpdef int vpoint_dof(object vpoints):
    """Degree of freedoms calculate from PMKS expressions."""
    # Joint with DOF 1.
    cdef int j1 = 0
    # Joint with DOF 2.
    cdef int j2 = 0
    # First link 'ground'.
    cdef set vlinks = {'ground'}

    cdef int link_count
    cdef VPoint vpoint
    for vpoint in vpoints:
        link_count = len(vpoint.links)
        if not link_count > 1:
            # If a point doesn't have two more links, it can not be call a 'joint'.
            continue
        vlinks.update(vpoint.links)
        if vpoint.type == VJoint.R:
            j1 += link_count - 1
        elif vpoint.type == VJoint.P:
            if link_count > 2:
                j1 += link_count - 2
            j1 += 1
        elif vpoint.type == VJoint.RP:
            if link_count > 2:
                j1 += link_count - 2
            j2 += 1
    return 3 * (len(vlinks) - 1) - (2 * j1) - j2


cdef inline int base_friend(int node, object vpoints):
    cdef int i
    cdef VPoint vpoint
    for i, vpoint in enumerate(vpoints):
        if not vpoints[node].links:
            continue
        if vpoints[node].links[0] in vpoint.links:
            return i


cdef inline double tuple_distance(tuple c1, tuple c2):
    """Calculate the distance between two tuple coordinates."""
    return distance(c1[0], c1[1], c2[0], c2[1])


cpdef tuple data_collecting(object exprs, dict mapping, object vpoints_):
    """Data collecting process.

    Input data:
    + exprs: [('PLAP', 'P0', 'L0', 'a0', 'P1', 'P2'), ...]
    + mapping: {0: 'P0', 1: 'P2', 2: 'P3', 3: 'P4', ...}
        + Specify position: {'P0': (10., 20.)}
        + Specify link length: {(a, b): 20.}
    + vpoints_: [VPoint0, VPoint1, VPoint2, ...]
    + pos: [(x0, y0), (x1, y1), (x2, y2), ...]
    
    vpoints will make a copy that we don't want to modified origin data.
    """
    cdef list vpoints = list(vpoints_)

    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix.
    cdef int node
    cdef str link
    cdef VPoint vpoint
    cdef dict vlinks = {}
    for node, vpoint in enumerate(vpoints):
        for link in vpoint.links:
            # Add as vlink.
            if link not in vlinks:
                vlinks[link] = [node]
            else:
                vlinks[link].append(node)

    # Replace the P joints and their friends with RP joint.
    # DOF must be same after properties changed.
    cdef int base
    cdef str link_
    cdef VPoint vpoint_
    cdef set links = set()
    for base in range(len(vpoints)):
        vpoint = vpoints[base]
        if vpoint.type != VJoint.P:
            continue
        for link in vpoint.links[1:]:
            links.clear()
            for node in vlinks[link]:
                vpoint_ = vpoints[node]
                if node == base or (vpoint_.type in {VJoint.P, VJoint.RP}):
                    continue
                links.update(vpoint_.links)
                vpoints[node] = VPoint.c_slider_joint(
                    [vpoint.links[0]] + [
                        link_ for link_ in vpoint_.links
                        if (link_ not in vpoint.links)
                    ],
                    VJoint.RP,
                    vpoint.angle,
                    vpoint_.cx,
                    vpoint_.cy
                )

    # Reverse mapping, exclude specified link length.
    cdef dict mapping_r = {}
    cdef dict length = {}
    cdef dict data_dict = {}

    cdef object k, v
    for k, v in mapping.items():
        if type(k) == int:
            mapping_r[v] = k
            if mapping[k] in mapping:
                data_dict[mapping[k]] = mapping[mapping[k]]
        elif type(k) == tuple:
            length[frozenset(k)] = v
            data_dict[k] = v

    cdef list pos = []
    for vpoint in vpoints:
        if vpoint.type == VJoint.R:
            pos.append(vpoint.c[0])
        else:
            pos.append(vpoint.c[1])

    cdef int i, bf
    cdef double angle
    # Add slider slot virtual coordinates.
    for i, vpoint in enumerate(vpoints):
        # PLPP dependents.
        if vpoint.type != VJoint.RP:
            continue
        bf = base_friend(i, vpoints)
        angle = radians(
            vpoint.angle -
            vpoint.slope_angle(vpoints[bf], 1, 0) +
            vpoint.slope_angle(vpoints[bf], 0, 0)
        )
        pos.append((vpoint.c[1][0] + cos(angle), vpoint.c[1][1] + sin(angle)))
        mapping_r[f'S{i}'] = len(pos) - 1

    # Add data to 'data_dict' and counting DOF.
    cdef int dof = 0
    cdef int target
    cdef tuple expr
    cdef frozenset pair
    for expr in exprs:
        node = mapping_r[expr[1]]
        target = mapping_r[expr[-1]]
        # Point 1: expr[1]
        if expr[1] not in data_dict:
            data_dict[expr[1]] = pos[mapping_r[expr[1]]]
        if expr[0] == 'PLAP':
            # Link 1: expr[2]
            pair = frozenset({node, target})
            if pair in length:
                data_dict[expr[2]] = length[pair]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            # Point 2: expr[4]
            if expr[4] not in data_dict:
                data_dict[expr[4]] = pos[mapping_r[expr[4]]]
            # Inputs
            dof += 1
        elif expr[0] == 'PLLP':
            # Link 1: expr[2]
            pair = frozenset({node, target})
            if pair in length:
                data_dict[expr[2]] = length[pair]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            # Link 2: expr[3]
            pair = frozenset({mapping_r[expr[4]], target})
            if pair in length:
                data_dict[expr[3]] = length[pair]
            else:
                data_dict[expr[3]] = tuple_distance(pos[mapping_r[expr[4]]], pos[target])
            # Point 2: expr[4]
            if expr[4] not in data_dict:
                data_dict[expr[4]] = pos[mapping_r[expr[4]]]
        elif expr[0] == 'PLPP':
            # Link 1: expr[2]
            pair = frozenset({node, target})
            if pair in length:
                data_dict[expr[2]] = length[pair]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            # Point 2:  expr[3]
            if expr[3] not in data_dict:
                data_dict[expr[3]] = pos[mapping_r[expr[3]]]
        elif expr[0] == 'PXY':
            # X: expr[2]
            if expr[2] in mapping:
                data_dict[expr[2]] = mapping[expr[2]]
            else:
                data_dict[expr[2]] = pos[target][0] - pos[node][0]
            # Y: expr[3]
            if expr[3] in mapping:
                data_dict[expr[3]] = mapping[expr[3]]
            else:
                data_dict[expr[3]] = pos[target][1] - pos[node][1]

    # Other grounded R joints.
    for i, vpoint in enumerate(vpoints):
        if vpoint.grounded() and vpoint.type == VJoint.R:
            data_dict[mapping[i]] = vpoint.c[0]

    return data_dict, dof


cpdef list expr_solving(
    object exprs,
    dict mapping,
    object vpoints,
    object angles = None
):
    """Solving function.
    
    + exprs: [('PLAP', 'P0', 'L0', 'a0', 'P1'), ...]
    TODO: Change another way to specify the lengths.
    + mapping: {0: 'P0', ..., (0, 1): 20.0, ...}
    + vpoints: [VPoint]
    + angles: [[a0]: a0, [a1]: a1, ...]
    """
    # Blank sequences.
    if angles is None:
        angles = []

    cdef dict data_dict
    cdef int dof_input
    data_dict, dof_input = data_collecting(exprs, mapping, vpoints)

    # Check input number.
    cdef int dof = vpoint_dof(vpoints)
    if dof_input > dof:
        raise ValueError(
            f"wrong number of input parameters: {dof_input} / {dof}"
        )

    # Reverse mapping, exclude specified link length.
    cdef object k, v
    cdef dict mapping_r = {v: k for k, v in mapping.items() if type(k) == int}

    # Check input pairs.
    cdef tuple expr
    for expr in exprs:
        if expr[0] == 'PLAP':
            if vpoints[mapping_r[expr[1]]].grounded() and vpoints[mapping_r[expr[-1]]].grounded():
                raise ValueError("wrong driver definition.")

    # Angles.
    cdef double a
    cdef int i
    for i, a in enumerate(angles):
        data_dict[f'a{i}'] = radians(a)

    # Solve
    if exprs:
        expr_parser(exprs, data_dict)

    cdef dict p_data_dict = {}
    cdef bint has_not_solved = False

    # Add coordinate of known points.
    for i in range(len(vpoints)):
        # {1: 'A'} vs {'A': (10., 20.)}
        if mapping[i] in data_dict:
            p_data_dict[i] = data_dict[mapping[i]]
        else:
            has_not_solved = True

    # Calling Sketch Solve kernel and try to get the result.
    cdef list solved_bfgs = []
    if has_not_solved:

        # Add specified link lengths.
        for k, v in data_dict.items():
            if type(k) == tuple:
                p_data_dict[k] = v

        # Solve
        try:
            solved_bfgs = vpoint_solving(vpoints, {}, p_data_dict)
        except ValueError:
            raise ValueError("result contains failure from sketch solve")

    # Format:
    # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
    # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
    cdef list solved_points = []
    cdef VPoint vpoint
    for i in range(len(vpoints)):
        vpoint = vpoints[i]
        if mapping[i] in data_dict:
            # These points has been solved.
            if isnan(data_dict[mapping[i]][0]):
                raise ValueError(f"result contains failure: Point{i}")
            if vpoint.type == VJoint.R:
                solved_points.append(data_dict[mapping[i]])
            else:
                solved_points.append((vpoint.c[0], data_dict[mapping[i]]))
        elif solved_bfgs:
            # These points solved by Sketch Solve.
            if vpoint.type == VJoint.R:
                solved_points.append(solved_bfgs[i])
            else:
                solved_points.append((solved_bfgs[i][0], solved_bfgs[i][1]))
        else:
            # No answer.
            if vpoint.type == VJoint.R:
                solved_points.append(vpoint.c[0])
            else:
                solved_points.append(vpoint.c)
    return solved_points
