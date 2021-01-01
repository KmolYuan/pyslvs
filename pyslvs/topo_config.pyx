# -*- coding: utf-8 -*-
# cython: language_level=3, cdivision=True, boundscheck=False, wraparound=False
# cython: initializedcheck=False, nonecheck=False

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libc.math cimport sin, cos, M_PI
from libcpp.map cimport map
from numpy import zeros, float64 as f64
from .expression cimport VJoint, VPoint, VLink
from .tinycadlib cimport (
    PXY, PPP, PLA, PLAP, PLLP, PLPP, PALP, P_LABEL, L_LABEL, I_LABEL,
    A_LABEL, S_LABEL
)

ctypedef vector[Sym] Inputs
ctypedef map[int, bint] Status


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


cdef class EStack:
    """Triangle solution stack, generated from [`t_config`](#t_config).
    It is pointless to call the constructor.
    """

    cdef void add_pxy(self, Sym c1, Sym v1, Sym v2, Sym t) nogil:
        cdef Expr e
        e.func = PXY
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_ppp(self, Sym c1, Sym c2, Sym c3, Sym t) nogil:
        cdef Expr e
        e.func = PPP
        e.c1 = c1
        e.c2 = c2
        e.c3 = c3
        e.target = t
        self.stack.push_back(e)

    cdef void add_pla(self, Sym c1, Sym v1, Sym v2, Sym t) nogil:
        cdef Expr e
        e.func = PLA
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_plap(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t) nogil:
        cdef Expr e
        e.func = PLAP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_pllp(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t) nogil:
        cdef Expr e
        e.func = PLLP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_plpp(self, Sym c1, Sym v1, Sym c2, Sym c3, Sym t, bint op) nogil:
        cdef Expr e
        e.func = PLPP
        e.c1 = c1
        e.v1 = v1
        e.c2 = c2
        e.c3 = c3
        e.target = t
        e.op = op
        self.stack.push_back(e)

    cdef void add_palp(self, Sym c1, Sym v1, Sym v2, Sym c2, Sym t, bint op) nogil:
        cdef Expr e
        e.func = PALP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.target = t
        e.op = op
        self.stack.push_back(e)

    cpdef list as_list(self):
        """Copy the dataset as list object."""
        stack = []
        cdef Expr expr
        for expr in self.stack:
            if expr.func == PXY:
                stack.append((
                    "PXY",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.target),
                ))
            elif expr.func == PPP:
                stack.append((
                    "PPP",
                    symbol_str(expr.c1),
                    symbol_str(expr.c2),
                    symbol_str(expr.c3),
                    symbol_str(expr.target),
                ))
            elif expr.func == PLA:
                stack.append((
                    "PLA",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.target),
                ))
            elif expr.func == PLAP:
                stack.append((
                    "PLAP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                    symbol_str(expr.target),
                ))
            elif expr.func == PLLP:
                stack.append((
                    "PLLP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                    symbol_str(expr.target),
                ))
            elif expr.func == PLPP:
                stack.append((
                    "PLPP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.c2),
                    symbol_str(expr.c3),
                    symbol_str(expr.target),
                ))
            elif expr.func == PALP:
                stack.append((
                    "PALP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                    symbol_str(expr.target),
                ))
        return stack

    def __repr__(self):
        return f"{type(self).__name__}({self.as_list()})"


cdef bint _is_all_lock(Status &status):
    """Test is all status done."""
    cdef bint n_status
    for _, n_status in status:
        if not n_status:
            return False
    return True


cdef bint _is_parallel(VPoint p1, VPoint p2, VPoint p3, VPoint p4):
    """Check parallel four bar loop."""
    return (abs(p1.distance(p3) - p2.distance(p4)) < 1e-12
            and abs(p2.distance(p3) - p1.distance(p4)) < 1e-12)


cdef bint _clockwise(double[:] c1, double[:] c2, double[:] c3):
    """Check the orientation of three points."""
    cdef double val = (c2[1] - c1[1]) * (c3[0] - c2[0]) - (c2[0] - c1[0]) * (c3[1] - c2[1])
    return val == 0 or val > 0


cdef (bint, int, int) _get_reliable_friends(
    int node,
    VPoint vpoint,
    object vlinks,
    Status &status
):
    """Return two "has been solved" nodes on the same link."""
    cdef int fa = -1
    cdef int fb = -1
    cdef int f
    for link in vpoint.links:
        if len(vlinks[link]) < 2:
            continue
        for f in vlinks[link]:
            if status[f] and f != node:
                if fa == -1:
                    fa = f
                elif fb == -1:
                    fb = f
                    return True, fa, fb
    return False, -1, -1


cdef int _get_intersection(
    int node1,
    int node2,
    object vpoints,
    object vlinks,
    Status &status
):
    """Get a configured node between two nodes."""
    cdef int node
    cdef VPoint vp
    for node, vp in enumerate(vpoints):
        if (set(vp.links) & set((<VPoint>vpoints[node1]).links)
            and set(vp.links) & set((<VPoint> vpoints[node2]).links)
            and status[node]
        ):
            return node
    return -1


cdef int _get_not_base_friend(VPoint vpoint, object vlinks, Status &status):
    """Get a configured node from other links."""
    if (vpoint.pin_grounded()
        or not vpoint.grounded()
        or vpoint.has_offset()
        or len(vpoint.links) < 2
    ):
        return -1
    cdef int f
    for f in vlinks[vpoint.links[1]]:
        if status[f]:
            return f
    return -1


cdef (bint, int, int) _get_base_friend(VPoint vpoint, object vlinks):
    """Get two configured nodes on the same links."""
    if len(vpoint.links) < 1:
        return False, -1, -1
    cdef int fa = -1
    cdef int f
    for f in vlinks[vpoint.links[0]]:
        if fa == -1:
            fa = f
        else:
            return True, fa, f
    return False, -1, -1


cdef int _get_input_base(int node, Inputs &inputs):
    """Get the base node for input pairs."""
    cdef int base, node_
    for base, node_ in inputs:
        if node == node_:
            return base
    return -1


cpdef EStack t_config(
    object vpoints_,
    object inputs_,
    object status_ = None
):
    """Generate the Triangle solution stack by mechanism expression `vpoints_`.

    The argument `inputs` is a list of input pairs.
    The argument `status` will track the configuration of each point,
    which is optional.
    """
    # For VPoint list:
    # + vpoints_: [vpoint0, vpoint1, ...]
    # + inputs: [(p0, p1), (p0, p2), ...]
    # + status: Dict[int, bint]
    # vpoints will make a copy that we don't want to modified itself
    cdef bint has_input = True
    cdef bint has_status = True
    if inputs_ is None or not inputs_:
        has_input = False
    if status_ is None:
        has_status = False
    cdef EStack exprs = EStack.__new__(EStack)
    if not vpoints_ or not has_input:
        return exprs

    vpoints = list(vpoints_)
    cdef Inputs inputs
    cdef Status status

    cdef bint ok
    cdef int node, base
    if has_input:
        for node, base in inputs_:
            inputs.push_back(Sym(node, base))
    if has_status:
        for node, ok in status_.items():
            status[node] = ok
    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix
    cdef VPoint vp1, vp2, vp3
    vlinks = {}
    for node, vp1 in enumerate(vpoints):
        status[node] = False
        if vp1.links:
            for link in vp1.links:
                # Connect on the ground and it is not a slider
                if link == VLink.FRAME and vp1.type == VJoint.R:
                    status[node] = True
                # Add as vlink
                if link not in vlinks:
                    vlinks[link] = {node}
                else:
                    vlinks[link].add(node)
        else:
            status[node] = True
    # Replace the P joints and their friends with RP joint
    # DOF must be same after properties changed
    for base, vp1 in enumerate(vpoints):
        if vp1.type != VJoint.P or not vp1.grounded():
            continue
        for link in vp1.links[1:]:
            for node in vlinks[link]:
                vp2 = vpoints[node]
                if node == base or vp2.type != VJoint.R:
                    continue
                vpoints[node] = VPoint.c_slider_joint([vp1.links[0]] + [
                    link_ for link_ in vp2.links
                    if link_ not in vp1.links
                ], VJoint.RP, vp1.angle, vp2.x, vp2.y)
    # Add positions parameters
    cdef double[:, :] pos = zeros((len(vpoints), 2), dtype=f64)
    for base, vp1 in enumerate(vpoints):
        node = 0 if vp1.type == VJoint.R else 1
        pos[base, 0] = vp1.c[node, 0]
        pos[base, 1] = vp1.c[node, 1]
    cdef int link_symbol = 0
    cdef int input_symbol = 0
    cdef int angle_symbol = 0
    # Input joints (R) that was connect with ground
    for base, node in inputs:
        if status[base]:
            exprs.add_pla(
                [P_LABEL, base],
                [L_LABEL, link_symbol],
                [I_LABEL, input_symbol],
                [P_LABEL, node]
            )
            status[node] = True
            link_symbol += 1
            input_symbol += 1
    # Now let we search around all the points, until we find the solutions
    # that we could
    input_targets = {node for _, node in inputs}
    node = 0
    cdef int skip_times = 0
    cdef int around = len(status)
    cdef int fa, fb, fc, fd
    cdef double angle
    cdef double[:] tmp = zeros(2, dtype=f64)
    # Friend iterator
    while not _is_all_lock(status):
        if status.find(node) == status.end():
            node = 0
            continue
        # Check and break the loop if it's re-scan again
        if skip_times >= around:
            break
        if status[node]:
            node += 1
            skip_times += 1
            continue
        vp1 = vpoints[node]
        if vp1.type == VJoint.R:
            # R joint
            # + Is input node?
            # + Normal revolute joint
            if node in input_targets:
                base = _get_input_base(node, inputs)
                if status[base]:
                    exprs.add_pla(
                        Sym(P_LABEL, base),
                        Sym(L_LABEL, link_symbol),
                        Sym(I_LABEL, input_symbol),
                        Sym(P_LABEL, node)
                    )
                    status[node] = True
                    link_symbol += 1
                    input_symbol += 1
                else:
                    skip_times += 1
            else:
                ok, fa, fb = _get_reliable_friends(node, vp1, vlinks, status)
                if not ok:
                    skip_times += 1
                else:
                    if not _clockwise(pos[fa], pos[node], pos[fb]):
                        fa, fb = fb, fa
                    vp2 = vpoints[fa]
                    vp3 = vpoints[fb]
                    if vp2.same_link(vp3) and not (vp2.grounded() and vp3.grounded()):
                        exprs.add_plap(
                            Sym(P_LABEL, fa),
                            Sym(L_LABEL, link_symbol),
                            Sym(A_LABEL, angle_symbol),
                            Sym(P_LABEL, fb),
                            Sym(P_LABEL, node)
                        )
                        link_symbol += 1
                        angle_symbol += 1
                    else:
                        # Decide when to solve parallel linkage
                        fc = _get_intersection(fa, fb, vpoints, vlinks, status)
                        if (fc != -1 and _is_parallel(
                            vpoints[fa], vpoints[fb], vpoints[fc], vpoints[node]
                        )):
                            exprs.add_ppp(
                                Sym(P_LABEL, fc),
                                Sym(P_LABEL, fa),
                                Sym(P_LABEL, fb),
                                Sym(P_LABEL, node)
                            )
                        else:
                            exprs.add_pllp(
                                Sym(P_LABEL, fa),
                                Sym(L_LABEL, link_symbol),
                                Sym(L_LABEL, link_symbol + 1),
                                Sym(P_LABEL, fb),
                                Sym(P_LABEL, node)
                            )
                            link_symbol += 2
                    status[node] = True
                    skip_times = 0
        elif vp1.type == VJoint.P:
            # Need to solve P joint itself here (only grounded)
            fa = _get_not_base_friend(vp1, vlinks, status)
            if fa == -1:
                skip_times += 1
            else:
                exprs.add_pxy(
                    Sym(P_LABEL, fa),
                    Sym(L_LABEL, link_symbol),
                    Sym(L_LABEL, link_symbol + 1),
                    Sym(P_LABEL, node)
                )
                status[node] = True
                link_symbol += 2
                # Solution for all friends
                for link in vp1.links[1:]:
                    for fb in vlinks[link]:
                        if status[fb]:
                            continue
                        exprs.add_pxy(
                            Sym(P_LABEL, node),
                            Sym(L_LABEL, link_symbol),
                            Sym(L_LABEL, link_symbol + 1),
                            Sym(P_LABEL, fb)
                        )
                        status[fb] = True
                        link_symbol += 2
                skip_times = 0
        elif vp1.type == VJoint.RP:
            # Copy as 'fc'
            fc = node
            # 'S' point
            tmp[:] = pos[node, :]
            angle = vp1.angle / 180 * M_PI
            tmp[0] += cos(angle)
            tmp[1] += sin(angle)
            fa = _get_not_base_friend(vp1, vlinks, status)
            ok, fb, fd = _get_base_friend(vp1, vlinks)
            if fa == -1 or not ok:
                skip_times += 1
            else:
                # Slot is not grounded
                if not vp1.grounded():
                    if not _clockwise(pos[fb], tmp, pos[fd]):
                        fb, fd = fd, fb
                    exprs.add_pllp(
                        Sym(P_LABEL, fb),
                        Sym(L_LABEL, link_symbol),
                        Sym(L_LABEL, link_symbol + 1),
                        Sym(P_LABEL, fd),
                        Sym(P_LABEL, node)
                    )
                    link_symbol += 2
                # PLPP
                # [PLLP]
                # Set 'S' (slider) point to define second point of slider
                # + A 'friend' from base link
                # + Get distance from me and friend
                # [PLPP]
                # Re-define coordinate of target point by self and 'S' point
                # + A 'friend' from other link
                # + Solve
                if not _clockwise(pos[fb], tmp, pos[fc]):
                    fb, fc = fc, fb
                exprs.add_pllp(
                    Sym(P_LABEL, fb),
                    Sym(L_LABEL, link_symbol),
                    Sym(L_LABEL, link_symbol + 1),
                    Sym(P_LABEL, fc),
                    Sym(S_LABEL, node)
                )
                # Two conditions
                exprs.add_plpp(
                    Sym(P_LABEL, fa),
                    Sym(L_LABEL, link_symbol + 2),
                    Sym(P_LABEL, node),
                    Sym(S_LABEL, node),
                    Sym(P_LABEL, node),
                    (pos[fa, 0] - pos[node, 0] > 0) != (vp1.angle > 90)
                )
                status[node] = True
                link_symbol += 3
                skip_times = 0
        node += 1
    if has_status:
        for node, ok in status:
            status_[node] = ok
    return exprs
