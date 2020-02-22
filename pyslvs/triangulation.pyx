# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions.

author: Yuan Chang
copyright: Copyright (C) 2016-2020
license: AGPL
email: pyslvs@gmail.com
"""

from typing import Sequence, Iterator
from libc.math cimport sin, cos, M_PI
from .expression cimport VJoint, VPoint, VLink


cdef inline str symbol_str(symbol p):
    """Pair to string."""
    if p.first == P_LABEL:
        return f"P{p.second}"
    elif p.first == L_LABEL:
        return f"L{p.second}"
    elif p.first == A_LABEL:
        return f"a{p.second}"
    elif p.first == S_LABEL:
        return f"S{p.second}"
    else:
        return ""


cdef class ExpressionStack:

    """Triangle solution stack, generated from
    [`vpoints_configure`](#vpoints_configure).
    It is pointless to call the constructor.
    """

    cdef void add_pla(self, symbol c1, symbol v1, symbol v2, symbol target):
        cdef Expression e
        e.func = PLA
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = target
        e.op = False
        self.stack.push_back(e)

    cdef void add_plap(self, symbol c1, symbol v1, symbol v2, symbol c2, symbol target):
        cdef Expression e
        e.func = PLAP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.c3 = target
        e.op = False
        self.stack.push_back(e)

    cdef void add_pllp(self, symbol c1, symbol v1, symbol v2, symbol c2, symbol target):
        cdef Expression e
        e.func = PLLP
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = c2
        e.c3 = target
        e.op = False
        self.stack.push_back(e)

    cdef void add_plpp(self, symbol c1, symbol v1, symbol c2, symbol c3, symbol target, bint op):
        cdef Expression e
        e.func = PLPP
        e.c1 = c1
        e.v1 = v1
        e.c2 = c2
        e.c3 = c3
        e.c4 = target
        e.op = op
        self.stack.push_back(e)

    cdef void add_pxy(self, symbol c1, symbol v1, symbol v2, symbol target):
        cdef Expression e
        e.func = PXY
        e.c1 = c1
        e.v1 = v1
        e.v2 = v2
        e.c2 = target
        e.op = False
        self.stack.push_back(e)

    cpdef list as_list(self):
        """Copy the dataset as list object."""
        stack = []
        cdef Expression expr
        for expr in self.stack:
            if expr.func == PLA:
                stack.append((
                    "PLAP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                ))
            elif expr.func == PLAP:
                stack.append((
                    "PLAP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                    symbol_str(expr.c3),
                ))
            elif expr.func == PLLP:
                stack.append((
                    "PLLP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                    symbol_str(expr.c3),
                ))
            elif expr.func == PLPP:
                stack.append((
                    "PLPP",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.c2),
                    symbol_str(expr.c3),
                    symbol_str(expr.c4),
                ))
            elif expr.func == PXY:
                stack.append((
                    "PXY",
                    symbol_str(expr.c1),
                    symbol_str(expr.v1),
                    symbol_str(expr.v2),
                    symbol_str(expr.c2),
                ))
        return stack

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.as_list()})"


cdef inline bint _is_all_lock(dict status):
    """Test is all status done."""
    cdef int node
    cdef bint n_status
    for node, n_status in status.items():
        if not n_status:
            return False
    return True


cdef inline bint _clockwise(tuple c1, tuple c2, tuple c3):
    """Check orientation of three points."""
    cdef double val = (c2[1] - c1[1]) * (c3[0] - c2[0]) - (c2[0] - c1[0]) * (c3[1] - c2[1])
    return val == 0 or val > 0


def _get_reliable_friend(
    node: int,
    vpoints: Sequence[VPoint],
    vlinks: dict,
    status: dict
) -> Iterator[int]:
    """Return a generator yield the vertices
        that "has been solved" on the same link.
    """
    cdef int friend
    for link in vpoints[node].links:
        if len(vlinks[link]) < 2:
            continue
        for friend in vlinks[link]:
            if status[friend] and friend != node:
                yield friend


def _get_not_base_friend(
    node: int,
    vpoints: Sequence[VPoint],
    vlinks: dict,
    status: dict
) -> Iterator[int]:
    """Get a friend from constrained nodes."""
    if len(vpoints[node].links) < 2:
        raise StopIteration
    cdef int friend
    for friend in vlinks[vpoints[node].links[1]]:
        if status[friend]:
            yield friend


def _get_base_friend(
    node: int,
    vpoints: Sequence[VPoint],
    vlinks: dict,
    status: dict
) -> Iterator[int]:
    """Get the constrained node of same links."""
    if len(vpoints[node].links) < 1:
        raise StopIteration
    cdef int friend
    for friend in vlinks[vpoints[node].links[0]]:
        yield friend


cdef inline int _get_input_base(int node, object inputs):
    """Get the base node for input pairs."""
    cdef int base, node_
    for base, node_ in inputs:
        if node == node_:
            return base
    return -1


cpdef ExpressionStack vpoints_configure(
    object vpoints_,
    object inputs,
    dict status = None
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

    if not vpoints_:
        return ExpressionStack.__new__(ExpressionStack)
    if not inputs:
        return ExpressionStack.__new__(ExpressionStack)

    vpoints = list(vpoints_)
    cdef int vpoints_count = len(vpoints)

    # First, we create a "VLinks" that can help us to
    # find a relationship just like adjacency matrix.
    cdef int node
    cdef VPoint vpoint
    vlinks = {}
    for node, vpoint in enumerate(vpoints):
        status[node] = False
        if vpoint.links:
            for link in vpoint.links:
                # Connect on the ground and it is not a slider.
                if link == VLink.FRAME and vpoint.type == VJoint.R:
                    status[node] = True
                # Add as vlink.
                if link not in vlinks:
                    vlinks[link] = {node}
                else:
                    vlinks[link].add(node)
        else:
            status[node] = True

    # Replace the P joints and their friends with RP joint.
    # DOF must be same after properties changed.
    cdef int base
    cdef VPoint vpoint_
    for base in range(vpoints_count):
        vpoint = vpoints[base]
        if vpoint.type != VJoint.P or not vpoint.grounded():
            continue
        for link in vpoint.links[1:]:
            links = set()
            for node in vlinks[link]:
                vpoint_ = vpoints[node]
                if node == base or vpoint_.type != VJoint.R:
                    continue
                links.update(vpoint_.links)
                vpoints[node] = VPoint.c_slider_joint(
                    [vpoint.links[0]] + [
                        link_ for link_ in vpoint_.links
                        if (link_ not in vpoint.links)
                    ],
                    VJoint.RP,
                    vpoint.angle,
                    vpoint_.x,
                    vpoint_.y
                )

    # Add positions parameters.
    pos = []
    for vpoint in vpoints:
        pos.append(vpoint.c[0 if vpoint.type == VJoint.R else 1])

    cdef ExpressionStack exprs = ExpressionStack.__new__(ExpressionStack)
    cdef int link_symbol = 0
    cdef int angle_symbol = 0

    # Input joints (R) that was connect with ground.
    for base, node in inputs:
        if status[base]:
            exprs.add_pla(
                [P_LABEL, base],
                [L_LABEL, link_symbol],
                [A_LABEL, angle_symbol],
                [P_LABEL, node]
            )
            status[node] = True
            link_symbol += 1
            angle_symbol += 1

    # Now let we search around all of points, until find the solutions that we could.
    input_targets = {node for base, node in inputs}
    node = 0
    cdef int skip_times = 0
    cdef int around = len(status)

    cdef int friend_a, friend_b, friend_c, friend_d
    cdef double tmp_x, tmp_y, angle
    # Friend iterator
    while not _is_all_lock(status):
        if node not in status:
            node = 0
            continue
        # Check and break the loop if it's re-scan again.
        if skip_times >= around:
            break
        if status[node]:
            node += 1
            skip_times += 1
            continue
        vpoint = vpoints[node]
        if vpoint.type == VJoint.R:
            # R joint
            # + Is input node?
            # + Normal revolute joint.
            if node in input_targets:
                base = _get_input_base(node, inputs)
                if status[base]:
                    exprs.add_pla(
                        [P_LABEL, base],
                        [L_LABEL, link_symbol],
                        [A_LABEL, angle_symbol],
                        [P_LABEL, node]
                    )
                    status[node] = True
                    link_symbol += 1
                    angle_symbol += 1
                else:
                    skip_times += 1
            else:
                fi = _get_reliable_friend(node, vpoints, vlinks, status)
                try:
                    friend_a = next(fi)
                    friend_b = next(fi)
                except StopIteration:
                    skip_times += 1
                else:
                    if not _clockwise(pos[friend_a], pos[node], pos[friend_b]):
                        friend_a, friend_b = friend_b, friend_a
                    exprs.add_pllp(
                        [P_LABEL, friend_a],
                        [L_LABEL, link_symbol],
                        [L_LABEL, link_symbol + 1],
                        [P_LABEL, friend_b],
                        [P_LABEL, node]
                    )
                    status[node] = True
                    link_symbol += 2
                    skip_times = 0

        elif vpoint.type == VJoint.P:
            # Need to solve P joint itself here (only grounded)
            fi = _get_not_base_friend(node, vpoints, vlinks, status)
            try:
                if vpoints[node].pin_grounded():
                    raise StopIteration
                if not vpoints[node].grounded():
                    raise StopIteration
                if vpoints[node].has_offset():
                    raise StopIteration
                friend_a = next(fi)
            except StopIteration:
                skip_times += 1
            else:
                exprs.add_pxy(
                    [P_LABEL, friend_a],
                    [L_LABEL, link_symbol],
                    [L_LABEL, link_symbol + 1],
                    [P_LABEL, node]
                )
                status[node] = True
                link_symbol += 2
                # Solution for all friends.
                for link in vpoints[node].links[1:]:
                    for friend_b in vlinks[link]:
                        if status[friend_b]:
                            continue
                        exprs.add_pxy(
                            [P_LABEL, node],
                            [L_LABEL, link_symbol],
                            [L_LABEL, link_symbol + 1],
                            [P_LABEL, friend_b]
                        )
                        status[friend_b] = True
                        link_symbol += 2
                skip_times = 0

        elif vpoint.type == VJoint.RP:
            # RP joint
            fi = _get_base_friend(node, vpoints, vlinks, status)
            # Copy as 'friend_c'.
            friend_c = node
            # 'S' point.
            tmp_x, tmp_y = pos[node]
            angle = vpoints[node].angle / 180 * M_PI
            tmp_x += cos(angle)
            tmp_y += sin(angle)
            try:
                if vpoints[node].pin_grounded():
                    raise StopIteration
                if not vpoints[node].grounded():
                    raise StopIteration
                if vpoints[node].has_offset():
                    raise StopIteration
                friend_a = next(_get_not_base_friend(node, vpoints, vlinks, status))
                friend_b = next(fi)
                # Slot is not grounded.
                if not vpoints[node].grounded():
                    friend_d = next(fi)
                    if not _clockwise(pos[friend_b], (tmp_x, tmp_y), pos[friend_d]):
                        friend_b, friend_d = friend_d, friend_b
                    exprs.add_pllp(
                        [P_LABEL, friend_b],
                        [L_LABEL, link_symbol],
                        [L_LABEL, link_symbol + 1],
                        [P_LABEL, friend_d],
                        [P_LABEL, node]
                    )
                    link_symbol += 2
            except StopIteration:
                skip_times += 1
            else:
                # PLPP
                # [PLLP]
                # Set 'S' (slider) point to define second point of slider.
                # + A 'friend' from base link.
                # + Get distance from me and friend.
                # [PLPP]
                # Re-define coordinate of target point by self and 'S' point.
                # + A 'friend' from other link.
                # + Solve.
                if not _clockwise(pos[friend_b], (tmp_x, tmp_y), pos[friend_c]):
                    friend_b, friend_c = friend_c, friend_b
                exprs.add_pllp(
                    [P_LABEL, friend_b],
                    [L_LABEL, link_symbol],
                    [L_LABEL, link_symbol + 1],
                    [P_LABEL, friend_c],
                    [S_LABEL, node]
                )
                # Two conditions.
                exprs.add_plpp(
                    [P_LABEL, friend_a],
                    [L_LABEL, link_symbol + 2],
                    [P_LABEL, node],
                    [S_LABEL, node],
                    [P_LABEL, node],
                    (pos[friend_a][0] - pos[node][0] > 0) != (vpoints[node].angle > 90)
                )
                status[node] = True
                link_symbol += 3
                skip_times = 0

        node += 1

    # exprs: [('PLAP', 'P0', 'L0', 'a0', 'P1', 'P2'), ...]
    return exprs
