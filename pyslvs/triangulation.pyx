# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from libc.math cimport sin, cos, M_PI
from numpy import zeros, float64 as np_float
from .expression cimport VJoint, VPoint, VLink


cdef str symbol_str(sym p):
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
    """Triangle solution stack, generated from
    [`t_config`](#t_config).
    It is pointless to call the constructor.
    """

    cdef void add_pla(self, sym c1, sym v1, sym v2, sym t) nogil:
        cdef Expr e
        e.func = PLA
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_plap(self, sym c1, sym v1, sym v2, sym c2, sym t) nogil:
        cdef Expr e
        e.func = PLAP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_pllp(self, sym c1, sym v1, sym v2, sym c2, sym t) nogil:
        cdef Expr e
        e.func = PLLP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cdef void add_plpp(self, sym c1, sym v1, sym c2, sym c3, sym t, bint op) nogil:
        cdef Expr e
        e.func = PLPP
        e.c1 = c1
        e.v1 = v1
        e.c2 = c2
        e.c3 = c3
        e.target = t
        e.op = op
        self.stack.push_back(e)

    cdef void add_pxy(self, sym c1, sym v1, sym v2, sym t) nogil:
        cdef Expr e
        e.func = PXY
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.target = t
        e.op = False
        self.stack.push_back(e)

    cpdef list as_list(self):
        """Copy the dataset as list object."""
        stack = []
        cdef Expr expr
        for expr in self.stack:
            if expr.func == PLA:
                stack.append((
                    "PLAP",
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
            elif expr.func == PXY:
                stack.append((
                    "PXY",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.target),
                ))
        return stack

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.as_list()})"


cdef bint _is_all_lock(object status):
    """Test is all status done."""
    cdef bint n_status
    for _, n_status in status.items():
        if not n_status:
            return False
    return True


@cython.boundscheck(False)
@cython.wraparound(False)
cdef bint _clockwise(double[:] c1, double[:] c2, double[:] c3):
    """Check orientation of three points."""
    cdef double val = (c2[1] - c1[1]) * (c3[0] - c2[0]) - (c2[0] - c1[0]) * (c3[1] - c2[1])
    return val == 0 or val > 0


cdef (bint, int, int) _get_reliable_friend(
    int node,
    VPoint vpoint,
    object vlinks,
    object status
):
    """Return a generator yield the vertices that "has been solved" on the
    same link.
    """
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
    return False, fa, fb


cdef int _get_not_base_friend(
    int node,
    VPoint vpoint,
    object vlinks,
    object status
):
    """Get a friend from constrained nodes."""
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


cdef (bint, int, int) _get_base_friend(
    int node,
    VPoint vpoint,
    object vlinks,
    object status
):
    """Get the constrained node of same links."""
    if len(vpoint.links) < 1:
        return False, -1, -1
    cdef int fa = -1
    cdef int fb = -1
    cdef int f
    for f in vlinks[vpoint.links[0]]:
        if fa == -1:
            fa = f
        elif fb == -1:
            fb = f
            return True, fa, fb
    return False, fa, fb


cdef int _get_input_base(int node, object inputs):
    """Get the base node for input pairs."""
    cdef int base, node_
    for base, node_ in inputs:
        if node == node_:
            return base
    return -1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef EStack t_config(
    object vpoints_,
    object inputs,
    object status = None
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
    if inputs is None:
        inputs = ()
    if status is None:
        status = {}
    cdef EStack exprs = EStack.__new__(EStack)
    if not vpoints_ or not inputs:
        return exprs
    vpoints = list(vpoints_)
    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix
    cdef int node
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
    cdef int base
    for base in range(len(vpoints)):
        vp1 = vpoints[base]
        if vp1.type != VJoint.P or not vp1.grounded():
            continue
        for link in vp1.links[1:]:
            links = set()
            for node in vlinks[link]:
                vp2 = vpoints[node]
                if node == base or vp2.type != VJoint.R:
                    continue
                links.update(vp2.links)
                vpoints[node] = VPoint.c_slider_joint([vp1.links[0]] + [
                    link_ for link_ in vp2.links
                    if link_ not in vp1.links
                ], VJoint.RP, vp1.angle, vp2.x, vp2.y)
    # Add positions parameters
    cdef double[:, :] pos = zeros((len(vpoints), 2), dtype=np_float)
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
    # Now let we search around all of points,
    # until find the solutions that we could
    input_targets = {node for _, node in inputs}
    node = 0
    cdef bint ok
    cdef int skip_times = 0
    cdef int around = len(status)
    cdef int fa, fb, fc, fd
    cdef double angle
    cdef double[:] tmp = zeros(2, dtype=np_float)
    # Friend iterator
    while not _is_all_lock(status):
        if node not in status:
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
                        sym(P_LABEL, base),
                        sym(L_LABEL, link_symbol),
                        sym(I_LABEL, input_symbol),
                        sym(P_LABEL, node)
                    )
                    status[node] = True
                    link_symbol += 1
                    input_symbol += 1
                else:
                    skip_times += 1
            else:
                ok, fa, fb = _get_reliable_friend(node, vp1, vlinks, status)
                if not ok:
                    skip_times += 1
                else:
                    if not _clockwise(pos[fa], pos[node], pos[fb]):
                        fa, fb = fb, fa
                    vp2 = vpoints[fa]
                    vp3 = vpoints[fb]
                    if vp2.same_link(vp3):
                        exprs.add_plap(
                            sym(P_LABEL, fa),
                            sym(L_LABEL, link_symbol),
                            sym(A_LABEL, angle_symbol),
                            sym(P_LABEL, fb),
                            sym(P_LABEL, node)
                        )
                        link_symbol += 1
                        angle_symbol += 1
                    else:
                        exprs.add_pllp(
                            sym(P_LABEL, fa),
                            sym(L_LABEL, link_symbol),
                            sym(L_LABEL, link_symbol + 1),
                            sym(P_LABEL, fb),
                            sym(P_LABEL, node)
                        )
                        link_symbol += 2
                    status[node] = True
                    skip_times = 0
        elif vp1.type == VJoint.P:
            # Need to solve P joint itself here (only grounded)
            fa = _get_not_base_friend(node, vp1, vlinks, status)
            if fa == -1:
                skip_times += 1
            else:
                exprs.add_pxy(
                    sym(P_LABEL, fa),
                    sym(L_LABEL, link_symbol),
                    sym(L_LABEL, link_symbol + 1),
                    sym(P_LABEL, node)
                )
                status[node] = True
                link_symbol += 2
                # Solution for all friends
                for link in vp1.links[1:]:
                    for fb in vlinks[link]:
                        if status[fb]:
                            continue
                        exprs.add_pxy(
                            sym(P_LABEL, node),
                            sym(L_LABEL, link_symbol),
                            sym(L_LABEL, link_symbol + 1),
                            sym(P_LABEL, fb)
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
            fa = _get_not_base_friend(node, vp1, vlinks, status)
            ok, fb, fd = _get_base_friend(node, vp1, vlinks, status)
            if fa == -1 or not ok:
                skip_times += 1
            else:
                # Slot is not grounded
                if not vp1.grounded():
                    if not _clockwise(pos[fb], tmp, pos[fd]):
                        fb, fd = fd, fb
                    exprs.add_pllp(
                        sym(P_LABEL, fb),
                        sym(L_LABEL, link_symbol),
                        sym(L_LABEL, link_symbol + 1),
                        sym(P_LABEL, fd),
                        sym(P_LABEL, node)
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
                    sym(P_LABEL, fb),
                    sym(L_LABEL, link_symbol),
                    sym(L_LABEL, link_symbol + 1),
                    sym(P_LABEL, fc),
                    sym(S_LABEL, node)
                )
                # Two conditions
                exprs.add_plpp(
                    sym(P_LABEL, fa),
                    sym(L_LABEL, link_symbol + 2),
                    sym(P_LABEL, node),
                    sym(S_LABEL, node),
                    sym(P_LABEL, node),
                    (pos[fa, 0] - pos[node, 0] > 0) != (vp1.angle > 90)
                )
                status[node] = True
                link_symbol += 3
                skip_times = 0
        node += 1
    return exprs
