# -*- coding: utf-8 -*-
# cython: language_level=3

"""Triangular expressions."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from typing import Sequence, Iterator
from libc.math cimport sin, cos
from tinycadlib cimport radians
from pmks cimport VPoint
from cpython cimport bool


cdef inline bool isAllLock(dict status, dict same = {}):
    """Test is all status done."""
    cdef int node
    cdef bool n_status
    for node, n_status in status.items():
        if (not n_status) and (node not in same):
            return False
    return True


cdef inline bool clockwise(tuple c1, tuple c2, tuple c3):
    """Check orientation of three points."""
    cdef double val = (c2[1] - c1[1])*(c3[0] - c2[0]) - (c2[0] - c1[0])*(c3[1] - c2[1])
    return ((val == 0) or (val > 0))


def _get_reliable_friend(
    node: int,
    vpoints: Sequence[VPoint],
    vlinks: dict,
    status: dict
) -> Iterator[int]:
    """Return a generator yield the nodes
        that "has been solved" on the same link.
    """
    cdef str link
    cdef int friend
    for link in vpoints[node].links:
        if len(vlinks[link]) < 2:
            continue
        for friend in vlinks[link]:
            if status[friend] and (friend != node):
                yield friend


def _get_notbase_friend(
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


cdef inline int get_input_base(int node, object inputs):
    """Get the base node for input pairs."""
    cdef int base, node_
    for base, node_ in inputs:
        if node == node_:
            return base
    return -1


cpdef list vpoints_configure(object vpoints_, object inputs = [], dict status = {}):
    """Auto configuration algorithm.
    
    For VPoint list.
    + vpoints_: [vpoint0, vpoint1, ...]
    + inputs: [(p0, p1), (p0, p2), ...]
    + status: Dict[int, bool]
    
    vpoints will make a copy that we don't want to modified itself.
    """
    if not vpoints_:
        return []
    
    cdef list vpoints = list(vpoints_)
    
    # First, we create a "VLinks" that can help us to
    #    find a releationship just like adjacency matrix.
    cdef int node
    cdef str link
    cdef VPoint vpoint
    cdef dict vlinks = {}
    for node, vpoint in enumerate(vpoints):
        status[node] = False
        if vpoint.links:
            for link in vpoint.links:
                # Connect on the ground and it is not a slider.
                if ('ground' == link) and (vpoint.type == 0):
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
    cdef str link_
    cdef VPoint vpoint_
    cdef set links
    for base in range(len(vpoints)):
        vpoint = vpoints[base]
        if (vpoint.type != 1) or (not vpoint.grounded()):
            continue
        for link in vpoint.links[1:]:
            links = set()
            for node in vlinks[link]:
                vpoint_ = vpoints[node]
                if (node == base) or (vpoint_.type != 0):
                    continue
                links.update(vpoint_.links)
                vpoints[node] = VPoint(
                    ",".join([vpoint.links[0]] + [
                        link_ for link_ in vpoint_.links
                        if (link_ not in vpoint.links)
                    ]),
                    2,
                    vpoint.angle,
                    vpoint_.colorSTR,
                    vpoint_.x,
                    vpoint_.y
                )
    
    # Add positions parameters.
    cdef list pos = []
    for vpoint in vpoints:
        pos.append(vpoint.c[0 if vpoint.type == 0 else 1])
    
    cdef list exprs = []
    cdef int link_symbol = 0
    cdef int angle_symbol = 0
    
    # Input joints (R) that was connect with ground.
    for base, node in inputs:
        if status[base]:
            exprs.append((
                'PLAP',
                f'P{base}',
                f'L{link_symbol}',
                f'a{angle_symbol}',
                f'P{node}',
            ))
            status[node] = True
            link_symbol += 1
            angle_symbol += 1
    
    # Now let we search around all of points, until find the solutions that we could.
    cdef set input_targets = {node for base, node in inputs}
    node = 0
    cdef int skip_times = 0
    cdef int around = len(status)
    
    cdef int type_num, friend_a, friend_b, friend_c, friend_d
    cdef double tmp_x, tmp_y, angle
    cdef bool reverse
    # Friend iterator.
    cdef object fi
    while not isAllLock(status):
        
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
        
        type_num = vpoints[node].type
        
        if type_num == 0:
            """R joint.
            
            + Is input node?
            + Normal revolute joint.
            """
            
            if node in input_targets:
                base = get_input_base(node, inputs)
                if status[base]:
                    exprs.append((
                        'PLAP',
                        f'P{base}',
                        f'L{link_symbol}',
                        f'a{angle_symbol}',
                        f'P{node}',
                    ))
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
                    if not clockwise(pos[friend_a], pos[node], pos[friend_b]):
                        friend_a, friend_b = friend_b, friend_a
                    exprs.append((
                        'PLLP',
                        f'P{friend_a}',
                        f'L{link_symbol}',
                        f'L{link_symbol + 1}',
                        f'P{friend_b}',
                        f'P{node}',
                    ))
                    status[node] = True
                    link_symbol += 2
                    skip_times = 0
        
        elif type_num == 1:
            """Need to solve P joint itself here. (only grounded)"""
            fi = _get_notbase_friend(node, vpoints, vlinks, status)
            try:
                if not vpoints[node].grounded():
                    raise StopIteration
                friend_a = next(fi)
            except StopIteration:
                skip_times += 1
            else:
                exprs.append((
                    'PXY',
                    f'P{friend_a}',
                    f'L{link_symbol}',
                    f'L{link_symbol + 1}',
                    f'P{node}',
                ))
                status[node] = True
                link_symbol += 2
                # Solution for all friends.
                for link in vpoints[node].links[1:]:
                    for friend_b in vlinks[link]:
                        if status[friend_b]:
                            continue
                        exprs.append((
                            'PXY',
                            f'P{node}',
                            f'L{link_symbol}',
                            f'L{link_symbol + 1}',
                            f'P{friend_b}',
                        ))
                        status[friend_b] = True
                        link_symbol += 2
                skip_times = 0
        
        elif type_num == 2:
            """RP joint."""
            fi = _get_base_friend(node, vpoints, vlinks, status)
            # Copy as 'friend_c'.
            friend_c = node
            # 'S' point.
            tmp_x, tmp_y = pos[node]
            angle = radians(vpoints[node].angle)
            tmp_x += cos(angle)
            tmp_y += sin(angle)
            try:
                friend_a = next(_get_notbase_friend(node, vpoints, vlinks, status))
                friend_b = next(fi)
                # Slot is not grounded.
                if not vpoints[node].grounded():
                    friend_d = next(fi)
                    if not clockwise(pos[friend_b], (tmp_x, tmp_y), pos[friend_d]):
                        friend_b, friend_d = friend_d, friend_b
                    exprs.append((
                        'PLLP',
                        f'P{friend_b}',
                        f'L{link_symbol}',
                        f'L{link_symbol + 1}',
                        f'P{friend_d}',
                        f'P{node}',
                    ))
                    link_symbol += 2
            except StopIteration:
                skip_times += 1
            else:
                """PLPP triangular.
                
                [PLLP]
                Set 'S' (slider) point to define second point of slider.
                + A 'friend' from base link.
                + Get distance from me and friend.
                
                [PLPP]
                Re-define coordinate of target point by self and 'S' point.
                + A 'friend' from other link.
                + Solving.
                """
                if not clockwise(pos[friend_b], (tmp_x, tmp_y), pos[friend_c]):
                    friend_b, friend_c = friend_c, friend_b
                exprs.append((
                    'PLLP',
                    f'P{friend_b}',
                    f'L{link_symbol}',
                    f'L{link_symbol + 1}',
                    f'P{friend_c}',
                    f'S{node}',
                ))
                # Two conditions.
                reverse = (pos[friend_a][0] - pos[node][0] > 0) != (vpoints[node].angle > 90)
                exprs.append((
                    'PLPP',
                    f'P{friend_a}',
                    f'L{link_symbol + 2}',
                    f'P{node}',
                    f'S{node}',
                    'T' if reverse else 'F',
                    f'P{node}',
                ))
                status[node] = True
                link_symbol += 3
                skip_times = 0
        
        node += 1
    
    # exprs: [('PLAP', 'P0', 'L0', 'a0', 'P1', 'P2'), ...]
    return exprs
