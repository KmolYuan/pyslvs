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
from .expression cimport Coord, VJoint, VPoint, VLink, distance, slope_angle
from .bfgs cimport SolverSystem


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


cpdef int vpoint_dof(object vpoints):
    """Return the DOF of the mechanism expression `vpoints`."""
    # Joint with DOF 1
    cdef int j1 = 0
    # Joint with DOF 2
    cdef int j2 = 0
    # First link is frame
    vlinks = {VLink.FRAME}
    cdef int link_count
    cdef VPoint vp
    for vp in vpoints:
        link_count = len(vp.links)
        if link_count < 2:
            # If a point doesn't have two more links, it can not be call a 'joint'.
            continue
        vlinks.update(vp.links)
        if vp.type in {VJoint.R, VJoint.P}:
            j1 += link_count - 1
        elif vp.type == VJoint.RP:
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


cdef bint preprocessing(EStack exprs, object vpoints, object angles,
                        map[Sym, CCoord] &joint_pos,
                        map[SwappablePair, double] &link_len,
                        map[Sym, double] &param):
    """Data extraction. Return true if input angle is matched DOF.

    Use "exprs", "vpoints", "angles" and "link_len" to generate solver
    required data. Please pre-allocate the "joint_pos" and "param".

    C++ objects
    + "exprs.stack" used for the position solution.
    + "joint_pos" used for joint position.
    + "link_len" used for custom link length except from calculation.
      (The link pairs must be exist in the solution.)
    + "param" used for store angles and link length.
    """
    vpoints = list(vpoints)
    cdef int vp_dof = vpoint_dof(vpoints)
    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix
    cdef int node, base
    cdef VPoint vp, vp2
    vlinks = {}
    for node, vp in enumerate(vpoints):
        for link in vp.links:
            # Add as vlink.
            if link not in vlinks:
                vlinks[link] = [node]
            else:
                vlinks[link].append(node)
    # Replace the P joints and their friends with RP joint
    # DOF must be same after properties changed
    cdef double x, y, angle
    for base in range(len(vpoints)):
        vp = vpoints[base]
        if vp.type != VJoint.P:
            continue
        for link in vp.links[1:]:
            for node in vlinks[link]:
                vp2 = vpoints[node]
                if node == base or vp2.is_slider():
                    continue
                x = vp2.c[0, 0]
                y = vp2.c[0, 1]
                vpoints[node] = VPoint.c_slider_joint([vp.links[0]] + [
                    link_ for link_ in vp2.links
                    if link_ not in vp.links
                ], VJoint.RP, vp.angle, x, y)
    # Assign joint positions
    for node, vp in enumerate(vpoints):
        base = 1 if vp.is_slider() else 0
        joint_pos[Sym(P_LABEL, node)] = CCoord(vp.c[base, 0], vp.c[base, 1])
    # Add slider slot virtual coordinates
    for node, vp in enumerate(vpoints):
        # PLPP dependencies
        if vp.type != VJoint.RP:
            continue
        base = base_friend(node, vpoints)
        angle = (vp.angle
                 - vp.slope_angle(vpoints[base], 1, 0)
                 + vp.slope_angle(vpoints[base], 0, 0)) / 180 * M_PI
        joint_pos[Sym(S_LABEL, node)] = CCoord(vp.c[1, 0] + cos(angle),
                                               vp.c[1, 1] + sin(angle))
    # Input angles
    for node, angle in enumerate(angles):
        param[Sym(I_LABEL, node)] = angle / 180 * M_PI
    # Scan again to check if the parameter exists
    # Especially link lengths and angles
    cdef int dof = 0
    cdef Expr e
    cdef SwappablePair pair1, pair2
    for e in exprs.stack:
        pair1 = SwappablePair(e.c1.second, e.target.second)
        pair2 = SwappablePair(e.c2.second, e.target.second)
        if e.func == PXY:
            param[e.v1] = joint_pos[e.target].x - joint_pos[e.c1].x
            param[e.v2] = joint_pos[e.target].y - joint_pos[e.c1].y
            continue
        if e.func in {PLLP, PALP}:
            param[e.v2] = (link_len[pair2]
                           if link_len.find(pair2) != link_len.end() else
                           distance(joint_pos[e.c2].x,
                                    joint_pos[e.c2].y,
                                    joint_pos[e.target].x,
                                    joint_pos[e.target].y))
        if e.func in {PLLP, PLA, PLAP, PLPP}:
            param[e.v1] = (link_len[pair1]
                           if link_len.find(pair1) != link_len.end() else
                           distance(joint_pos[e.c1].x,
                                    joint_pos[e.c1].y,
                                    joint_pos[e.target].x,
                                    joint_pos[e.target].y))
            if e.func in {PLA, PLAP}:
                if e.v2.first == I_LABEL:
                    dof += 1
                else:  # A_LABEL
                    param[e.v2] = slope_angle(joint_pos[e.c1].x,
                                              joint_pos[e.c1].y,
                                              joint_pos[e.target].x,
                                              joint_pos[e.target].y)
                    if e.func == PLAP:
                        param[e.v2] -= slope_angle(joint_pos[e.c1].x,
                                                   joint_pos[e.c1].y,
                                                   joint_pos[e.c2].x,
                                                   joint_pos[e.c2].y)
                if (
                    (<VPoint>vpoints[e.c1.second]).grounded()
                    and (<VPoint>vpoints[e.target.second]).grounded()
                ):
                    raise ValueError("wrong driver definition")
    return len(angles) == dof <= vp_dof


cpdef list expr_solving(
    EStack exprs,
    object vpoints,
    object angles = None
):
    """Solver function of Triangular method and BFGS method, for mechanism
    expression `vpoints`.

    The triangle expression stack `expr` is generated from
    [`t_config`](#t_config).

    Solver function will not handle slider input pairs in argument `angles`,
    which is only support revolute joints. In another way, the slider input
    pairs can be set by [`VPoint.disable_offset()`](#vpointdisable_offset)
    method.
    """
    # Blank sequences
    if angles is None:
        angles = []
    cdef map[Sym, CCoord] joint_pos
    cdef map[SwappablePair, double] link_len
    cdef map[Sym, double] param
    if not preprocessing(exprs, vpoints, angles, joint_pos, link_len, param):
        raise ValueError("wrong number of input parameters")
    # Check coverage
    status = {i: (<VPoint> vp).grounded() for i, vp in enumerate(vpoints)}
    cdef Expr e
    for e in exprs.stack:
        status[e.target.second] = True
    cdef bint bfgs_mode = not all(status.values())
    # Solve
    cdef ExprSolver solver = ExprSolver(exprs.stack, joint_pos, link_len, param)
    if not solver.solve():
        raise ValueError("solve failed")
    # Use BFGS mode
    cdef pair[Sym, CCoord] jp
    if bfgs_mode:
        data_dict = {}
        for jp in solver.joint_pos:
            if jp.first.first == P_LABEL and status[jp.first.second]:
                data_dict[jp.first.second] = Coord.__new__(Coord, jp.second.x,
                                                           jp.second.y)
        bfgs_rt = SolverSystem(vpoints, {}, data_dict).solve()
    rt = []
    cdef int i
    cdef CCoord c
    cdef VPoint vp
    for i, vp in enumerate(vpoints):
        if bfgs_mode:
            if vp.is_slider():
                rt.append((bfgs_rt[i][0], bfgs_rt[i][1]))
            else:
                rt.append(bfgs_rt[i])
        else:
            c = solver.joint_pos[Sym(P_LABEL, i)]
            if vp.is_slider():
                rt.append(((vp.c[0, 0], vp.c[0, 1]), (c.x, c.y)))
            else:
                rt.append((c.x, c.y))
    return rt

cdef (bint, map[Sym, CCoord]) quick_solve(
    vector[Expr] stack,
    map[Sym, CCoord] joint_pos,
    map[SwappablePair, double] link_len,
    map[Sym, double] param
) nogil:
    """Quick solving function.

    !! Cython can not expose C external function to another pyx.
    """
    cdef ExprSolver s = ExprSolver(stack, joint_pos, link_len, param)
    cdef bint ok = s.solve()
    return ok, s.joint_pos
