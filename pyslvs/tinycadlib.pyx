# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Tiny CAD library of PMKS symbolic and position analysis.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from libc.math cimport M_PI, sin, cos
from .expression cimport Coord, VJoint, VPoint, VLink
from .bfgs cimport SolverSystem


cdef inline double radians(double degree) nogil:
    """Deg to rad."""
    return degree / 180 * M_PI


cdef str symbol_str(Sym p):
    """The match pattern of the symbols."""
    if p.first == P_LABEL:
        return f"P{p.second}"
    elif p.first == L_LABEL:
        return f"L{p.second}"
    elif p.first == I_LABEL:
        return f"I{p.second}"
    elif p.first == A_LABEL:
        return f"A{p.second}"
    elif p.first == S_LABEL:
        return f"S{p.second}"
    else:
        return ""


def pxy(Coord c1, double x, double y):
    """The PXY function requires one point and offset values, get the
    position of second point.

    In the following picture, `c1` correspond to "A", `d0` correspond to "X",
    `d1` correspond to "Y", `return` correspond to "B", the sign of value are
    correspond to coordinate system.

    ![pxy](img/pxy.png)
    """
    cdef CCoord c = cpxy(CCoord(c1.x, c1.y), x, y)
    return Coord.__new__(Coord, c.x, c.y)


def ppp(Coord c1, Coord c2, Coord c3):
    """The PPP function is used to solve parallel linkage.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `c3` correspond to "C", `return` correspond to "D".

    ![ppp](img/ppp.png)
    """
    cdef CCoord c = cppp(CCoord(c1.x, c1.y), CCoord(c2.x, c2.y),
                         CCoord(c3.x, c3.y))
    return Coord.__new__(Coord, c.x, c.y)


def plap(
    Coord c1,
    double d0,
    double a0,
    Coord c2 = None,
    bint inverse = False
):
    """The PLAP function requires two points, one distance and one angle,
    obtained the position of third point. The unit of `a0` is degree.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `d0` correspond to "L0", `a0` correspond to "beta", `return` correspond
    to "C".
    If `c2` is not given, "alpha" will be set to zero.

    ![plap](img/plap.png)

    Set `inverse` option to `True` can make `a0` value as negative.
    """
    if c2 is None:
        c2 = c1
    cdef CCoord c = cplap(CCoord(c1.x, c1.y), d0, a0, CCoord(c2.x, c2.y), inverse)
    return Coord.__new__(Coord, c.x, c.y)


def pllp(
    Coord c1,
    double d0,
    double d1,
    Coord c2,
    bint inverse = False
):
    """The PLLP function requires two points and two distances, obtained the
    position of third point.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `d0` correspond to "L0", `d1` correspond to "L1", `return` correspond to
    "C".

    ![pllp](img/pllp.png)

    Set `inverse` option to `True` can make the result upside down.
    """
    cdef CCoord c = cpllp(CCoord(c1.x, c1.y), d0, d1, CCoord(c2.x, c2.y), inverse)
    return Coord.__new__(Coord, c.x, c.y)


def plpp(
    Coord c1,
    double d0,
    Coord c2,
    Coord c3,
    bint inverse = False
):
    """The PLPP function requires three points and one distance, obtained the
    position of fourth point.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `c3` correspond to "C", `d0` correspond to "L0", `return` correspond to "D".

    ![plpp](img/plpp.png)

    Set `inverse` option to `True` can make the result to the another side
    between `c1` and line `c2` `c3`.
    """
    cdef CCoord c = cplpp(CCoord(c1.x, c1.y), d0, CCoord(c2.x, c2.y),
                          CCoord(c3.x, c3.y), inverse)
    return Coord.__new__(Coord, c.x, c.y)


def palp(
    Coord c1,
    double a0,
    double d0,
    Coord c2,
    bint inverse = False
):
    """The PALP function requires two points, one angle and one distance,
    obtained the position of fourth point.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `d0` correspond to "L0", `a0` correspond to "alpha", `return` correspond
    to "C".

    ![palp](img/palp.png)

    Set `inverse` option to `True` can make the result upside down.
    """
    cdef CCoord c = cpalp(CCoord(c1.x, c1.y), a0, d0, CCoord(c2.x, c2.y), inverse)
    return Coord.__new__(Coord, c.x, c.y)


cpdef void expr_parser(EStack exprs, dict data_dict):
    """Solve and update information of the triangle expression `exprs` to 
    `data_dict`.
    The argument `exprs` can be obtained by
    [`t_config`](#t_config) and [`EStack.as_list()`](#estackas_list) method.

    This function is already included in [`expr_solving`](#expr_solving),
    not recommended for direct use.
    """
    # Update data
    # + exprs: [("PLAP", "P0", "L0", "a0", "P1", "P2"), ..."]
    # + data_dict: {'a0':0., 'L1':10., 'A':(30., 40.), ...}
    cdef Sym target
    cdef Coord coord
    cdef Expr expr
    for expr in exprs.stack:
        if expr.func == PXY:
            coord = pxy(
                data_dict[symbol_str(expr.c1)],
                data_dict[symbol_str(expr.v1)],
                data_dict[symbol_str(expr.v2)]
            )
        elif expr.func == PPP:
            coord = ppp(
                data_dict[symbol_str(expr.c1)],
                data_dict[symbol_str(expr.c2)],
                data_dict[symbol_str(expr.c3)]
            )
        elif expr.func in {PLA, PLAP}:
            if expr.func == PLA:
                coord = plap(
                    data_dict[symbol_str(expr.c1)],
                    data_dict[symbol_str(expr.v1)],
                    data_dict[symbol_str(expr.v2)]
                )
            else:
                coord = plap(
                    data_dict[symbol_str(expr.c1)],
                    data_dict[symbol_str(expr.v1)],
                    data_dict[symbol_str(expr.v2)],
                    data_dict[symbol_str(expr.c2)],
                    expr.op
                )
        elif expr.func == PLLP:
            coord = pllp(
                data_dict[symbol_str(expr.c1)],
                data_dict[symbol_str(expr.v1)],
                data_dict[symbol_str(expr.v2)],
                data_dict[symbol_str(expr.c2)],
                expr.op
            )
        elif expr.func == PLPP:
            coord = plpp(
                data_dict[symbol_str(expr.c1)],
                data_dict[symbol_str(expr.v1)],
                data_dict[symbol_str(expr.c2)],
                data_dict[symbol_str(expr.c3)],
                expr.op
            )
        elif expr.func == PALP:
            coord = palp(
                data_dict[symbol_str(expr.c1)],
                data_dict[symbol_str(expr.v1)],
                data_dict[symbol_str(expr.v2)],
                data_dict[symbol_str(expr.c2)],
                expr.op
            )
        else:
            raise ValueError("unsupported function")
        data_dict[symbol_str(expr.target)] = coord


cpdef int vpoint_dof(object vpoints):
    """Return the DOF of the mechanism expression `vpoints`."""
    # Joint with DOF 1
    cdef int j1 = 0
    # Joint with DOF 2
    cdef int j2 = 0
    # First link is frame
    vlinks = {VLink.FRAME}
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
    return 3 * (len(vlinks) - 1) - 2 * j1 - j2


cdef inline int base_friend(int node, object vpoints):
    cdef int i
    cdef VPoint vpoint
    for i, vpoint in enumerate(vpoints):
        if not vpoints[node].links:
            continue
        if vpoints[node].links[0] in vpoint.links:
            return i


cpdef int data_collecting(dict data_dict, EStack exprs, dict mapping,
                          object vpoints_):
    """Data transform function of Triangular method.
    The triangle expression stack `expr` is generated from
    [`t_config`](#t_config).
    The information data `mapping` map the symbols to the indicator of 
    `vpoints_`.

    This function is already included in [`expr_solving`](#expr_solving),
    not recommended for direct use.
    """
    vpoints = list(vpoints_)
    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix
    cdef int node
    cdef VPoint vpoint
    vlinks = {}
    for node, vpoint in enumerate(vpoints):
        for link in vpoint.links:
            # Add as vlink.
            if link not in vlinks:
                vlinks[link] = [node]
            else:
                vlinks[link].append(node)
    # Replace the P joints and their friends with RP joint
    # DOF must be same after properties changed
    cdef int base
    cdef double x, y
    cdef VPoint vpoint_
    links = set()
    for base in range(len(vpoints)):
        vpoint = vpoints[base]
        if vpoint.type != VJoint.P:
            continue
        for link in vpoint.links[1:]:
            links.clear()
            for node in vlinks[link]:
                vpoint_ = vpoints[node]
                if node == base or vpoint_.type in {VJoint.P, VJoint.RP}:
                    continue
                links.update(vpoint_.links)
                x = vpoint_.c[0, 0]
                y = vpoint_.c[0, 1]
                vpoints[node] = VPoint.c_slider_joint(
                    [vpoint.links[0]] + [
                        link_ for link_ in vpoint_.links
                        if (link_ not in vpoint.links)
                    ],
                    VJoint.RP,
                    vpoint.angle,
                    x,
                    y
                )
    # Reverse mapping, exclude specified link length
    mapping_r = {}
    length = {}
    for k, v in mapping.items():
        if type(k) is int:
            mapping_r[v] = k
            if v in mapping:
                x, y = mapping[v]
                data_dict[v] = Coord.__new__(Coord, x, y)
        elif type(k) is tuple:
            length[frozenset(k)] = v
    pos = []
    for vpoint in vpoints:
        if vpoint.type == VJoint.R:
            pos.append(Coord.__new__(Coord, vpoint.c[0, 0], vpoint.c[0, 1]))
        else:
            pos.append(Coord.__new__(Coord, vpoint.c[1, 0], vpoint.c[1, 1]))
    cdef int i, bf
    cdef double angle
    # Add slider slot virtual coordinates
    for i, vpoint in enumerate(vpoints):
        # PLPP dependencies
        if vpoint.type != VJoint.RP:
            continue
        bf = base_friend(i, vpoints)
        angle = radians(vpoint.angle
                        - vpoint.slope_angle(vpoints[bf], 1, 0)
                        + vpoint.slope_angle(vpoints[bf], 0, 0))
        pos.append(Coord.__new__(Coord, vpoint.c[1, 0] + cos(angle),
                                 vpoint.c[1, 1] + sin(angle)))
        mapping_r[symbol_str(Sym(S_LABEL, i))] = len(pos) - 1
    # Add data to 'data_dict' and counting DOF
    cdef int dof = 0
    cdef int target
    cdef Expr expr
    cdef Coord coord1, coord2
    for expr in exprs.stack:
        node = mapping_r[symbol_str(expr.c1)]
        base = mapping_r[symbol_str(expr.c2)]
        target = mapping_r[symbol_str(expr.target)]
        # Point 1
        if symbol_str(expr.c1) not in data_dict:
            data_dict[symbol_str(expr.c1)] = pos[node]
        # For all expressions
        if expr.func == PXY:
            coord1 = pos[node]
            coord2 = pos[target]
            # X
            if symbol_str(expr.v1) in mapping:
                x, y = mapping[symbol_str(expr.v1)]
                data_dict[symbol_str(expr.v1)] = Coord.__new__(Coord, x, y)
            else:
                data_dict[symbol_str(expr.v1)] = coord2.x - coord1.x
            # Y
            if symbol_str(expr.v2) in mapping:
                x, y = mapping[symbol_str(expr.v2)]
                data_dict[symbol_str(expr.v2)] = Coord.__new__(Coord, x, y)
            else:
                data_dict[symbol_str(expr.v2)] = coord2.y - coord1.y
        elif expr.func == PPP:
            # Point 2
            if symbol_str(expr.c2) not in data_dict:
                data_dict[symbol_str(expr.c2)] = pos[base]
            # Point 3
            if symbol_str(expr.c3) not in data_dict:
                data_dict[symbol_str(expr.c3)] = pos[mapping_r[symbol_str(expr.c3)]]
        elif expr.func in {PLA, PLAP}:
            # Link 1
            pair = frozenset({node, target})
            coord1 = pos[node]
            coord2 = pos[target]
            if pair in length:
                data_dict[symbol_str(expr.v1)] = length[pair]
            else:
                data_dict[symbol_str(expr.v1)] = coord1.distance(coord2)
            # Angle 1
            if expr.v2.first == I_LABEL:
                # Inputs
                dof += 1
            else:
                # Links
                if expr.func == PLA:
                    data_dict[symbol_str(expr.v2)] = coord1.slope_angle(coord2)
                else:
                    data_dict[symbol_str(expr.v2)] = (coord1.slope_angle(coord2)
                        - coord1.slope_angle(pos[base]))
            # Point 2
            if expr.func == PLAP and symbol_str(expr.c2) not in data_dict:
                data_dict[symbol_str(expr.c2)] = pos[base]
        elif expr.func == PLLP:
            # Link 1
            pair = frozenset({node, target})
            if pair in length:
                data_dict[symbol_str(expr.v1)] = length[pair]
            else:
                coord1 = pos[node]
                coord2 = pos[target]
                data_dict[symbol_str(expr.v1)] = coord1.distance(coord2)
            # Link 2
            pair = frozenset({base, target})
            if pair in length:
                data_dict[symbol_str(expr.v2)] = length[pair]
            else:
                coord1 = pos[base]
                coord2 = pos[target]
                data_dict[symbol_str(expr.v2)] = coord1.distance(coord2)
            # Point 2
            if symbol_str(expr.c2) not in data_dict:
                data_dict[symbol_str(expr.c2)] = pos[base]
        elif expr.func == PLPP:
            # Link 1
            pair = frozenset({node, target})
            if pair in length:
                data_dict[symbol_str(expr.v1)] = length[pair]
            else:
                coord1 = pos[node]
                coord2 = pos[target]
                data_dict[symbol_str(expr.v1)] = coord1.distance(coord2)
            # Point 2
            if symbol_str(expr.c2) not in data_dict:
                data_dict[symbol_str(expr.c2)] = pos[base]
            # Point 3 is not needed
        elif expr.func == PALP:
            # Angle 1
            vpoint = vpoints[node]
            coord1 = pos[node]
            coord2 = pos[base]
            data_dict[symbol_str(expr.v1)] = (vpoint.angle
                                              - coord1.slope_angle(coord2))
            # Link 1
            pair = frozenset({base, target})
            if pair in length:
                data_dict[symbol_str(expr.v2)] = length[pair]
            else:
                coord1 = pos[base]
                coord2 = pos[target]
                data_dict[symbol_str(expr.v2)] = coord1.distance(coord2)
            # Point 2
            data_dict[symbol_str(expr.c2)] = pos[base]
    # Other grounded R joints
    for i, vpoint in enumerate(vpoints):
        if vpoint.grounded() and vpoint.type == VJoint.R:
            x = vpoint.c[0, 0]
            y = vpoint.c[0, 1]
            data_dict[mapping[i]] = Coord.__new__(Coord, x, y)
    return dof


cpdef list expr_solving(
    EStack exprs,
    object vpoints,
    object angles = None
):
    """Solver function of Triangular method and BFGS method, for mechanism 
    expression `vpoints`.

    The triangle expression stack `expr` is generated from
    [`t_config`](#t_config).

    The information data `mapping` map the symbols to the indicator of 
    `vpoints`,
    additionally has a same format as argument `data_dict` in [SolverSystem].

    Solver function will not handle slider input pairs in argument `angles`,
    which is only support revolute joints. In another way, the slider input 
    pairs
    can be set by [`VPoint.disable_offset()`](#vpointdisable_offset) method.
    """
    # Blank sequences
    if angles is None:
        angles = []
    data_dict = {}
    mapping = {n: f'P{n}' for n in range(len(vpoints))}
    cdef int dof_input = data_collecting(data_dict, exprs, mapping, vpoints)
    # Check input number
    cdef int dof = vpoint_dof(vpoints)
    if dof_input > dof:
        raise ValueError(
            f"wrong number of input parameters: {dof_input} / {dof}"
        )
    # Reverse mapping, exclude specified link length
    mapping_r = {v: k for k, v in mapping.items() if type(k) is int}
    # Check input pairs
    cdef int target
    cdef Expr expr
    for expr in exprs.stack:
        if expr.func in {PLA, PLAP}:
            target = mapping_r[symbol_str(expr.target)]
            if (
                vpoints[mapping_r[symbol_str(expr.c1)]].grounded()
                and vpoints[target].grounded()
            ):
                raise ValueError("wrong driver definition.")
    # Angles
    cdef double a
    cdef int i
    for i, a in enumerate(angles):
        data_dict[symbol_str(Sym(I_LABEL, i))] = radians(a)
    # Solve
    if not exprs.stack.empty():
        expr_parser(exprs, data_dict)
    p_data_dict = {}
    cdef bint has_not_solved = False
    # Add coordinate of known points
    for i in range(len(vpoints)):
        # {1: 'A'} vs {'A': (10., 20.)}
        if mapping[i] in data_dict:
            p_data_dict[i] = data_dict[mapping[i]]
        else:
            has_not_solved = True
    # Calling Sketch Solve kernel and try to get the result
    if has_not_solved:
        # Add specified link lengths
        for k, v in mapping.items():
            if type(k) is tuple:
                p_data_dict[k] = v
        # Solve
        try:
            solved_bfgs = SolverSystem(vpoints, {}, p_data_dict).solve()
        except ValueError:
            raise ValueError("result contains failure from sketch solve")
    # Format:
    # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
    # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
    solved_points = []
    cdef VPoint vpoint
    cdef Coord coord
    for i in range(len(vpoints)):
        vpoint = vpoints[i]
        if mapping[i] in data_dict:
            # These points has been solved
            coord = data_dict[mapping[i]]
            if coord.is_nan():
                raise ValueError(f"result contains failure: Point{i}")
            if vpoint.type == VJoint.R:
                solved_points.append((coord.x, coord.y))
            else:
                solved_points.append((
                    (vpoint.c[0, 0], vpoint.c[0, 1]),
                    (coord.x, coord.y)
                ))
        elif solved_bfgs is not None:
            # These points solved by Sketch Solve
            if vpoint.type == VJoint.R:
                solved_points.append(solved_bfgs[i])
            else:
                solved_points.append((solved_bfgs[i][0], solved_bfgs[i][1]))
        else:
            # No answer
            if vpoint.type == VJoint.R:
                solved_points.append(vpoint.c[0, :])
            else:
                solved_points.append(vpoint.c)
    return solved_points
