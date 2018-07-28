# -*- coding: utf-8 -*-
# cython: language_level=3

"""Tiny CAD library of PMKS simbolic and position analysis."""

# __author__ = "Yuan Chang"
# __copyright__ = "Copyright (C) 2016-2018"
# __license__ = "AGPL"
# __email__ = "pyslvs@gmail.com"

from cpython cimport bool
from libc.math cimport (
    M_PI,
    sqrt,
    sin,
    cos,
    atan2,
    hypot,
)
#Not a number
cdef double nan = float('nan')
from numpy import isnan
from pmks cimport VPoint
from bfgs cimport vpoint_solving


cdef inline double distance(double x1, double y1, double x2, double y2):
    """Distance of two cartesian coordinates."""
    return hypot(x2 - x1, y2 - y1)


cdef class Coordinate:
    
    """A class to store the coordinate."""
    
    def __cinit__(self, double x, double y):
        self.x = x
        self.y = y
    
    cpdef double distance(self, Coordinate p):
        """Distance."""
        return distance(self.x, self.y, p.x, p.y)
    
    cpdef bool is_nan(self):
        """Test this coordinate is a error-occured answer."""
        return bool(isnan(self.x))
    
    def __repr__(self):
        """Debug printing."""
        return f"Coordinate({self.x:.02f}, {self.y:.02f})"


cpdef tuple PLAP(
    Coordinate A,
    double L0,
    double a0,
    Coordinate B = None,
    bool inverse = False
):
    """Point on circle by angle."""
    cdef double a1 = atan2(B.y - A.y, B.x - A.x) if B else 0
    if inverse:
        return (A.x + L0 * cos(a1 - a0), A.y + L0 * sin(a1 - a0))
    else:
        return (A.x + L0 * cos(a1 + a0), A.y + L0 * sin(a1 + a0))


cpdef tuple PLLP(
    Coordinate A,
    double L0,
    double L1,
    Coordinate B,
    bool inverse = False
):
    """Two intersection points of two circles."""
    cdef double dx = B.x - A.x
    cdef double dy = B.y - A.y
    cdef double d = A.distance(B)
    
    #No solutions, the circles are separate.
    if d > L0 + L1:
        return (nan, nan)
    
    #No solutions because one circle is contained within the other.
    if d < abs(L0 - L1):
        return (nan, nan)
    
    #Circles are coincident and there are an infinite number of solutions.
    if (d == 0) and (L0 == L1):
        return (nan, nan)
    cdef double a = (L0 * L0 - L1 * L1 + d * d) / (2 * d)
    cdef double h = sqrt(L0 * L0 - a * a)
    cdef double xm = A.x + a * dx / d
    cdef double ym = A.y + a * dy / d
    
    if inverse:
        return (xm + h * dy / d, ym - h * dx / d)
    else:
        return (xm - h * dy / d, ym + h * dx / d)


cpdef tuple PLPP(
    Coordinate A,
    double L0,
    Coordinate B,
    Coordinate C,
    bool inverse = False
):
    """Two intersection points of a line and a circle."""
    cdef double line_mag = B.distance(C)
    cdef double dx = C.x - B.x
    cdef double dy = C.y - B.y
    cdef double u = ((A.x - B.x) * dx + (A.y - B.y) * dy) / (line_mag * line_mag)
    cdef Coordinate I = Coordinate(B.x + u * dx, B.y + u * dy)
    
    #Test distance between point A and intersection.
    cdef double d = A.distance(I)
    if d > L0:
        #No intersection.
        return (nan, nan)
    elif d == L0:
        #One intersection point.
        return (I.x, I.y)
    
    #Two intersection points.
    d = sqrt(L0 * L0 - d * d) / line_mag
    if inverse:
        return (I.x - dx * d, I.y - dy * d)
    else:
        return (I.x + dx * d, I.y + dy * d)


cpdef tuple PXY(Coordinate A, double x, double y):
    """Using relative cartesian coordinate to get solution."""
    return (A.x + x, A.y + y)


cdef inline bool legal_crank(Coordinate A, Coordinate B, Coordinate C, Coordinate D):
    """
    verify the fourbar is satisfied the Gruebler's Equation, s + g <= p + q
        C - D
        |   |
        A   B
    """
    cdef double driver = A.distance(C)
    cdef double follower = B.distance(D)
    cdef double ground = A.distance(B)
    cdef double connector = C.distance(D)
    return (
        (driver + connector <= ground + follower) or
        (driver + ground <= connector + follower)
    )


cdef inline str strbetween(str s, str front, str back):
    """Get the string that is inside of parenthesis."""
    return s[s.find(front) + 1:s.find(back)]


cdef inline str strbefore(str s, str front):
    """Get the string that is front of parenthesis."""
    return s[:s.find(front)]


cpdef void expr_parser(object exprs, dict data_dict):
    """Update data.
    
    + exprs: "PLAP[P0,L0,a0,P1](P2);PLLP[P2,L1,L2,P1](P3);..."
    + data_dict: {'a0':0., 'L1':10., 'A':(30., 40.), ...}
    """
    cdef object expr
    cdef str fun, p
    cdef list params
    cdef int params_count
    cdef double x, y
    for expr in exprs:
        #If the mechanism has no any solution.
        if not expr:
            return
        fun = expr[0]
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
        
        #We should unpack as C++'s way.
        x = float('nan')
        y = x
        if fun == 'PLAP':
            if params_count == 3:
                x, y = PLAP(params[0], params[1], params[2])
            elif params_count == 4:
                x, y = PLAP(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = PLAP(params[0], params[1], params[2], params[3], params[4])
        elif fun == 'PLLP':
            if params_count == 4:
                x, y = PLLP(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = PLLP(params[0], params[1], params[2], params[3], params[4])
        elif fun == 'PLPP':
            if params_count == 4:
                x, y = PLPP(params[0], params[1], params[2], params[3])
            elif params_count == 5:
                x, y = PLPP(params[0], params[1], params[2], params[3], params[4])
        elif fun == 'PXY':
            if params_count == 3:
                x, y = PXY(params[0], params[1], params[2])
        data_dict[expr[-1]] = (x, y)


cpdef str expr_join(object exprs):
    """Use to append a list of symbols into a string."""
    return ';'.join([
        f"{expr[0]}[{','.join(expr[1:-1]), expr[-1])}]({expr[-1]})" for expr in exprs
    ])


cpdef tuple expr_parse(str exprs):
    """Parse expression as tuple."""
    exprs = exprs.replace(" ", '')
    cdef list tmp_list = []
    cdef list params
    cdef str expr, p
    for expr in exprs.split(';'):
        if not expr:
            return ()
        params = []
        params.append(strbefore(expr, '['))
        params.extend(strbetween(expr, '[', ']').split(','))
        params.append(strbetween(expr, '(', ')'))
        tmp_list.append(tuple(params))
    return tuple(tmp_list)


cpdef int vpoint_dof(object vpoints):
    """Degree of freedoms calculate from PMKS expressions."""
    cdef int j1 = 0
    cdef int j2 = 0
    cdef set vlinks = {'ground'}
    
    cdef int link_count
    cdef VPoint vpoint
    for vpoint in vpoints:
        link_count = len(vpoint.links)
        if not link_count > 1:
            #If a point doesn't have two more links, it can not be call a 'joint'.
            continue
        vlinks.update(vpoint.links)
        if vpoint.type == 0:
            j1 += link_count - 1
        elif vpoint.type == 1:
            if link_count > 2:
                j1 += link_count - 2
            j1 += 1
        elif vpoint.type == 2:
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
    """Input data:
    
    + exprs: [('PLAP', 'P0', 'L0', 'a0', 'P1', 'P2'), ...]
    + mapping: {0: 'P0', 1: 'P2', 2: 'P3', 3: 'P4', ...}
        + Specify link length: mapping['L0'] = 20.0
    + vpoints_: [VPoint0, VPoint1, VPoint2, ...]
    + pos: [(x0, y0), (x1, y1), (x2, y2), ...]
    
    vpoints will make a copy that we don't want to modified origin data.
    """
    cdef list vpoints = list(vpoints_)
    
    #First, we create a "VLinks" that can help us to
    #find a releationship just like adjacency matrix.
    cdef int node
    cdef str link
    cdef VPoint vpoint
    cdef dict vlinks = {}
    for node, vpoint in enumerate(vpoints):
        for link in vpoint.links:
            #Add as vlink.
            if link not in vlinks:
                vlinks[link] = [node]
            else:
                vlinks[link].append(node)
    
    #Replace the P joints and their friends with RP joint.
    #DOF must be same after properties changed.
    cdef int base
    cdef str link_
    cdef VPoint vpoint_
    cdef set links = set()
    for base in range(len(vpoints)):
        vpoint = vpoints[base]
        if vpoint.type != 1:
            continue
        for link in vpoint.links[1:]:
            links.clear()
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
                    vpoint_.cx,
                    vpoint_.cy
                )
    
    #Reverse mapping, exclude specified link length.
    cdef object k, v
    cdef dict mapping_r = {v: k for k, v in mapping.items() if (type(k) == int)}
    
    cdef list pos = []
    for vpoint in vpoints:
        if vpoint.type == 0:
            pos.append(vpoint.c[0])
        else:
            pos.append(vpoint.c[1])
    
    cdef int i, bf
    cdef double angle
    #Add slider slot virtual coordinates.
    for i, vpoint in enumerate(vpoints):
        #PLPP dependents.
        if vpoint.type != 2:
            continue
        bf = base_friend(i, vpoints)
        angle = (
            vpoint.angle -
            vpoint.slope_angle(vpoints[bf], 1, 0) +
            vpoint.slope_angle(vpoints[bf], 0, 0)
        ) / 180 * M_PI
        pos.append((vpoint.c[1][0] + cos(angle), vpoint.c[1][1] + sin(angle)))
        mapping_r[f'S{i}'] = len(pos) - 1
    
    cdef int dof = 0
    cdef dict data_dict = {}
    cdef int target
    cdef tuple expr
    """Add data to 'data_dict'.
    
    + Add 'L' (link) parameters.
    + Counting DOF and targets.
    """
    for expr in exprs:
        node = mapping_r[expr[1]]
        target = mapping_r[expr[-1]]
        #Point 1: expr[1]
        if expr[1] not in data_dict:
            data_dict[expr[1]] = pos[mapping_r[expr[1]]]
        if expr[0] == 'PLAP':
            #Link 1: expr[2]
            if expr[2] in mapping:
                data_dict[expr[2]] = mapping[expr[2]]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            #Point 2: expr[4]
            if expr[4] not in data_dict:
                data_dict[expr[4]] = pos[mapping_r[expr[4]]]
            #Inputs
            dof += 1
        elif expr[0] == 'PLLP':
            #Link 1: expr[2]
            if expr[2] in mapping:
                data_dict[expr[2]] = mapping[expr[2]]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            #Link 2: expr[3]
            if expr[3] in mapping:
                data_dict[expr[3]] = mapping[expr[3]]
            else:
                data_dict[expr[3]] = tuple_distance(pos[mapping_r[expr[4]]], pos[target])
            #Point 2: expr[4]
            if expr[4] not in data_dict:
                data_dict[expr[4]] = pos[mapping_r[expr[4]]]
        elif expr[0] == 'PLPP':
            #Link 1: expr[2]
            if expr[2] in mapping:
                data_dict[expr[2]] = mapping[expr[2]]
            else:
                data_dict[expr[2]] = tuple_distance(pos[node], pos[target])
            #Point 2:  expr[3]
            if expr[3] not in data_dict:
                data_dict[expr[3]] = pos[mapping_r[expr[3]]]
        elif expr[0] == 'PXY':
            #X: expr[2]
            if expr[2] in mapping:
                data_dict[expr[2]] = mapping[expr[2]]
            else:
                data_dict[expr[2]] = pos[target][0] - pos[node][0]
            #Y: expr[3]
            if expr[3] in mapping:
                data_dict[expr[3]] = mapping[expr[3]]
            else:
                data_dict[expr[3]] = pos[target][1] - pos[node][1]
    #Other grounded R joints.
    for i, vpoint in enumerate(vpoints):
        if vpoint.grounded() and (vpoint.type == 0):
            data_dict[mapping[i]] = vpoint.c[0]
    return data_dict, dof


cpdef list expr_solving(
    object exprs,
    dict mapping,
    object vpoints,
    object angles = []
):
    """Solving function.
    
    + angles: [[a0]: a0, [a1]: a1, ...]
    """
    cdef dict data_dict
    cdef int dof_input
    data_dict, dof_input = data_collecting(exprs, mapping, vpoints)
    
    #Check input number.
    cdef int dof = vpoint_dof(vpoints)
    if dof_input > dof:
        raise Exception(
            f"wrong number of input parameters: {dof_input} / {dof}"
        )
    
    #Angles.
    cdef double a
    cdef int i
    for i, a in enumerate(angles):
        data_dict[f'a{i}'] = a / 180 * M_PI
    
    #Solve
    expr_parser(exprs, data_dict)
    
    cdef dict p_data_dict = {}
    cdef bool has_target = False
    for i in range(len(vpoints)):
        if mapping[i] in data_dict:
            p_data_dict[i] = data_dict[mapping[i]]
        else:
            has_target = True
            break
    cdef list solved_bfgs = []
    if exprs and has_target:
        try:
            solved_bfgs = vpoint_solving(vpoints, [], p_data_dict)
        except Exception as e:
            raise Exception("result contains failure: Sketch Solve") from e
    
    """Format:
    
    + R joint
        [p0]: (p0_x, p0_y)
        [p1]: (p1_x, p1_y)
    + P or RP joint
        [p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))
    """
    cdef list solved_points = []
    for i in range(len(vpoints)):
        if mapping[i] not in data_dict:
            #These points solved by Sketch Solve.
            if solved_bfgs:
                if vpoints[i].type == 0:
                    solved_points.append(solved_bfgs[i])
                else:
                    solved_points.append((vpoints[i].c[0], solved_bfgs[i]))
            else:
                if vpoints[i].type == 0:
                    solved_points.append(vpoints[i].c[0])
                else:
                    solved_points.append(vpoints[i].c)
        else:
            if isnan(data_dict[mapping[i]][0]):
                raise Exception(f"result contains failure: Point{i}")
            if vpoints[i].type == 0:
                solved_points.append(data_dict[mapping[i]])
            else:
                solved_points.append((vpoints[i].c[0], data_dict[mapping[i]]))
    return solved_points
