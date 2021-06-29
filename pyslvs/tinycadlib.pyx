# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Tiny CAD library of PMKS symbolic and position analysis.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libc.math cimport M_PI, NAN, sin, cos
from .expression cimport Coord, VJoint, VPoint, VLink, distance, slope_angle
from .bfgs cimport SolverSystem
from numpy import zeros, array
from numpy.random import uniform


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


cdef bint preprocessing(EStack exprs, object vpoints, object inputs,
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
                    # Input angles
                    angle = inputs[pair1.first, pair1.second]
                    param[e.v2] = angle / 180 * M_PI
                    dof += 1
                else:
                    # A_LABEL
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
    return dof <= len(inputs) <= vp_dof


cpdef list expr_solving(
    EStack exprs,
    object vpoints,
    object inputs = None
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
    if inputs is None:
        inputs = {}
    cdef map[Sym, CCoord] joint_pos
    cdef map[SwappablePair, double] link_len
    cdef map[Sym, double] param
    if not preprocessing(exprs, vpoints, inputs, joint_pos, link_len, param):
        raise ValueError("wrong number of input parameters")
    # Check coverage
    status = {i: False for i in range(len(vpoints))}
    cdef Expr e
    for e in exprs.stack:
        status[e.c1.second] = status[e.target.second] = True
        if e.func in {PLAP, PLLP, PLPP, PALP}:
            status[e.c2.second] = True
            if e.func == PLPP:
                status[e.c3.second] = True
    cdef bint bfgs_mode = not all(status.values())
    # Solve
    cdef ExprSolver solver = ExprSolver(exprs.stack, joint_pos, param)
    if not solver.solve():
        raise ValueError("solve failed")
    # Use BFGS mode
    cdef pair[Sym, CCoord] jp
    if bfgs_mode:
        data_dict = {}
        for jp in solver.joint_pos:
            if (
                jp.first.first == P_LABEL
                and status[jp.first.second]
            ):
                data_dict[jp.first.second] = Coord.__new__(Coord, jp.second.x,
                                                           jp.second.y)
        return SolverSystem(vpoints, inputs, data_dict).solve()
    rt = []
    cdef int i
    cdef CCoord c
    cdef VPoint vp
    for i, vp in enumerate(vpoints):
        c = solver.joint_pos[Sym(P_LABEL, i)]
        if vp.is_slider():
            rt.append(((vp.c[0, 0], vp.c[0, 1]), (c.x, c.y)))
        else:
            rt.append((c.x, c.y))
    return rt


cdef (bint, map[Sym, CCoord]) quick_solve(
    vector[Expr] stack,
    map[Sym, CCoord] joint_pos,
    map[Sym, double] param
) nogil:
    """Quick solving function.

    !! Cython can not expose C external function to another pyx.
    """
    cdef ExprSolver s = ExprSolver(stack, joint_pos, param)
    cdef bint ok = s.solve()
    return ok, s.joint_pos


def uniform_four_bar(double ml, int n):
    """Generate n four bar mechanisms from maximum lengths.

    These mechanisms have coupling points.
    Normalized parameters are $[L_0, L_2, L_3, L_4, \\alpha]$.

    ![pxy](img/uniform_four_bar.png)
    """
    return array(_uniform_four_bar(ml, n))


cdef double[:, :] _uniform_four_bar(double ml, int n):
    """Uniform four-bar implementation."""
    cdef double[:, :] d = uniform(1., ml, (n, 5))
    cdef int i
    for i in range(n):
        while not d[i, 0] + d[i, 2] > 1 + d[i, 1]:
            for j in range(4):
                d[i, j] = uniform(1., ml)
        d[i, 4] = uniform(1., 2 * M_PI)
    return d


def uniform_path(double[:, :] v, int n):
    """Generate path with four-bar dimensions.

    Normalized parameters are $[L_0, L_2, L_3, L_4, \\alpha]$.
    """
    return array(c_uniform_path(v, n))


cdef (vector[Expr], map[Sym, CCoord], map[Sym, double]) ueb(double[:] v) nogil:
    """Uniform expression builder."""
    cdef vector[Expr] stack
    stack.push_back(Expr(False, PLA, Sym(L_LABEL, 1), Sym(I_LABEL, 0),
                         Sym(P_LABEL, 0), Sym(), Sym(), Sym(P_LABEL, 2)))
    stack.push_back(Expr(False, PLLP, Sym(L_LABEL, 2), Sym(L_LABEL, 3),
                         Sym(P_LABEL, 2), Sym(P_LABEL, 1), Sym(),
                         Sym(P_LABEL, 3)))
    stack.push_back(Expr(False, PLAP, Sym(L_LABEL, 4), Sym(A_LABEL, 0),
                         Sym(P_LABEL, 2), Sym(P_LABEL, 3), Sym(),
                         Sym(P_LABEL, 4)))
    cdef map[Sym, CCoord] joint_pos
    joint_pos[Sym(P_LABEL, 0)] = CCoord(0, 0)
    cdef map[Sym, double] param
    param[Sym(L_LABEL, 1)] = 1.
    joint_pos[Sym(P_LABEL, 1)] = CCoord(v[0], 0)
    cdef int i
    for i in range(1, 4):
        param[Sym(L_LABEL, i + 1)] = v[i]
    param[Sym(A_LABEL, 0)] = v[4]
    # Assign driver angle
    # param[Sym(I_LABEL, 0)] = a
    return stack, joint_pos, param


cdef double[:, :, :] c_uniform_path(double[:, :] v, int n) nogil:
    """Uniform path implementation."""
    cdef double[:, :, :] p
    with gil:
        p = zeros((v.shape[0], n, 2))
    cdef vector[Expr] stack
    cdef map[Sym, CCoord] joint_pos
    cdef map[Sym, double] param
    cdef bint ok
    cdef int i, j
    cdef double a
    cdef CCoord c
    cdef map[Sym, CCoord] ans
    for i in range(len(v)):
        stack, joint_pos, param = ueb(v[i])
        j = 0
        a = 0
        while a < 2 * M_PI:
            param[Sym(I_LABEL, 0)] = a
            ok, ans = quick_solve(stack, joint_pos, param)
            if ok:
                c = ans[Sym(P_LABEL, 4)]
                p[i, j, 0] = c.x
                p[i, j, 1] = c.y
            else:
                p[i, j, 0] = NAN
                p[i, j, 1] = NAN
            a += 2 * M_PI / n
            j += 1
    return p


cpdef object uniform_expr(double[:] v):
    """Turn the uniform link length into expression."""
    cdef vector[Expr] stack
    cdef map[Sym, CCoord] joint_pos
    cdef map[Sym, double] param
    stack, joint_pos, param = ueb(v)
    cdef vector[CCoord] coords = vector[CCoord](3, CCoord(0., 0.))
    cdef bint ok = False
    cdef double a = 0
    cdef map[Sym, CCoord] ans
    while a < 2 * M_PI:
        param[Sym(I_LABEL, 0)] = a
        ok, ans = quick_solve(stack, joint_pos, param)
        for i in range(2, 5):
            coords[i - 2] = ans[Sym(P_LABEL, i)]
        if ok:
            break
        a += 2 * M_PI / 30
    return [
        VPoint.c_r_joint([VLink.FRAME, 'L1'], 0., 0.),
        VPoint.c_r_joint([VLink.FRAME, 'L3'], v[0], 0.),
        VPoint.c_r_joint(['L1', 'L2'], coords[0].x, coords[0].y),
        VPoint.c_r_joint(['L2', 'L3'], coords[1].x, coords[1].y),
        VPoint.c_r_joint(['L2'], coords[2].x, coords[2].y),
    ]
