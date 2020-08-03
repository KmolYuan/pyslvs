cimport cython
from collections import OrderedDict
from numpy import zeros, array
from libc.math cimport HUGE_VAL, M_PI
from .topo_config cimport (
    t_config, symbol_str, I_LABEL, A_LABEL, Expr, PXY, PPP, PLA, PLAP, PLLP,
    PLPP, PALP, EStack,
)
from .bfgs cimport SolverSystem
from .tinycadlib cimport radians, pxy, ppp, plap, pllp, plpp, palp
from .expression cimport VJoint, VPoint
from .metaheuristics.utility cimport ObjFunc


@cython.final
@cython.boundscheck(False)
@cython.wraparound(False)
cdef class FMatch(ObjFunc):
    """This class is used to verified kinematics of the linkage mechanism.

    A fast matching method that adds mapping angles to variables.
    """
    cdef bint bfgs_mode, shape_only, ordered
    cdef int target_count, input_count, l_base
    cdef list vpoints, mapping_list
    cdef dict placement, target, mapping, mapping_r, data_dict
    cdef object inputs
    cdef double[:] polar_angles
    cdef EStack exprs

    def __cinit__(self, dict mech):
        # mech = {
        #     'expression': List[VPoint],
        #     'input': OrderedDict([((b0, d0), [start, end]), ...]),
        #     'placement': {pt: (x, y, r)},
        #     'target': {pt: [(x0, y0), (x1, y1), ...]},
        #     'same': {pt: match_to_pt},
        #     # Bounds of base link length
        #     'upper': float,
        #     'lower': float,
        #     'shape_only': bool,
        # }
        placement = mech.get('placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")
        target = mech.get('target', {})
        if len(target) == 0:
            raise ValueError("no target joint")
        check_set = {len(t) for t in target.values()}
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()
        # Change the target paths into memory view
        self.target = {}
        same = mech.get('same', {})
        self.shape_only = mech.get('shape_only', False)
        cdef int i, j
        cdef Coord[:] path
        for i in target:
            path = array([Coord(x, y) for x, y in target[i]], dtype=object)
            for j in range(i):
                if j in same:
                    i -= 1
            if self.shape_only:
                _norm(path, 1)
            self.target[i] = path
        # Expressions
        self.vpoints = list(mech.get('expression', []))
        self.inputs = OrderedDict(mech.get('input', {}))
        status = {}
        self.exprs = t_config(self.vpoints, tuple(self.inputs.keys()), status)
        self.bfgs_mode = not all(status.values())
        # Data mapping
        self.mapping = {i: f"P{i}" for i in range(len(self.vpoints))}
        self.mapping_r = {v: k for k, v in self.mapping.items()}
        self.mapping_list = []
        # Bounds
        upper = []
        lower = []
        if len(upper) != len(lower):
            raise ValueError("upper and lower should be in the same size")
        # Position
        cdef double x, y, r
        for i in sorted(placement):
            x, y, r = placement[i]
            upper.append(x + r)
            upper.append(y + r)
            lower.append(x - r)
            lower.append(y - r)
            self.mapping_list.append(i)
        cdef int p_base = len(upper)
        # Length of links
        link_upper = float(mech.get('upper', 100))
        link_lower = float(mech.get('lower', 0))
        i = 0
        cdef Expr expr
        for expr in self.exprs.stack:
            upper.append(link_upper)
            lower.append(link_lower)
            sym = frozenset({expr.c1.second, expr.target.second})
            self.mapping[sym] = None
            self.mapping_list.append(sym)
            if expr.func in {PLA, PLAP}:
                if expr.v2.first == A_LABEL:
                    upper.append(2 * M_PI)
                    lower.append(0)
                    sym = f"A{i}"
                    self.mapping[sym] = None
                    self.mapping_list.append(sym)
                    i += 1
            elif expr.func == PLLP:
                upper.append(link_upper)
                lower.append(link_lower)
                sym = frozenset({expr.c2.second, expr.target.second})
                self.mapping[sym] = None
                self.mapping_list.append(sym)
        self.l_base = len(upper)
        # Input nodes
        self.input_count = len(self.inputs)
        for start, end in self.inputs.values():
            upper.append(radians(start))
            lower.append(radians(end))
        # Angle rage
        upper[self.l_base:] *= self.target_count
        lower[self.l_base:] *= self.target_count
        self.upper = array(upper, dtype=np_float)
        self.lower = array(lower, dtype=np_float)
        # Swap upper and lower bound if reversed
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]
        # Allocate memory
        self.polar_angles = zeros(self.l_base - p_base, dtype=np_float)
        # Result list
        self.data_dict = {}

    cpdef double[:] get_upper(self):
        """Return upper bound."""
        return self.upper

    cpdef double[:] get_lower(self):
        """Return lower bound."""
        return self.lower

    cpdef bint is_two_kernel(self):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        return self.bfgs_mode

    cdef inline double get_len(self, str expr1, str expr2):
        return self.mapping[frozenset({self.mapping_r[expr1], self.mapping_r[
            expr2]})]

    cdef bint solve(self, double[:] input_list):
        self.data_dict.clear()
        cdef int i
        cdef VPoint vpoint
        cdef Coord coord1, coord2, coord3
        for i, vpoint in enumerate(self.vpoints):
            coord1 = Coord.__new__(Coord, vpoint.c[0, 0], vpoint.c[0, 1])
            if vpoint.type == VJoint.R:
                self.data_dict[self.mapping[i]] = coord1
                self.data_dict[i, -1] = coord1
            else:
                coord2 = Coord.__new__(Coord, vpoint.c[1, 0], vpoint.c[1, 1])
                self.data_dict[self.mapping[i]] = coord2
                self.data_dict[i, -1] = coord1
                self.data_dict[i, -2] = coord2
        # Solve
        i = 0
        cdef int a = 0
        cdef int t, params_count
        cdef double length, angle
        cdef Expr expr
        for expr in self.exprs.stack:
            target = symbol_str(expr.target)
            if expr.func == PXY:
                coord1 = self.data_dict[symbol_str(expr.c1)]
                coord3 = pxy(
                    coord1,
                    vpoint.c[0, 0] - coord1.x,
                    vpoint.c[0, 1] - coord1.y
                )
            elif expr.func == PPP:
                coord3 = ppp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.data_dict[symbol_str(expr.c2)],
                    self.data_dict[symbol_str(expr.c3)],
                )
            elif expr.func in {PLA, PLAP}:
                coord1 = self.data_dict[symbol_str(expr.c1)]
                if expr.func == PLAP:
                    coord2 = self.data_dict[symbol_str(expr.c2)]
                length = self.get_len(symbol_str(expr.c1), target)
                if expr.v2.first == I_LABEL:
                    angle = input_list[i]
                    i += 1
                else:
                    angle = self.polar_angles[a]
                    a += 1
                if expr.func == PLA:
                    coord3 = plap(coord1, length, angle)
                else:
                    coord3 = plap(coord1, length, angle, coord2, expr.op)
            elif expr.func == PLLP:
                coord3 = pllp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.get_len(symbol_str(expr.c1), target),
                    self.get_len(symbol_str(expr.c2), target),
                    self.data_dict[symbol_str(expr.c2)],
                    expr.op
                )
            elif expr.func == PLPP:
                coord3 = plpp(
                    self.data_dict[symbol_str(expr.c1)],
                    self.get_len(symbol_str(expr.c1), target),
                    self.data_dict[symbol_str(expr.c2)],
                    self.data_dict[symbol_str(expr.c3)],
                    expr.op
                )
            elif expr.func == PALP:
                angle = self.polar_angles[a]
                a += 1
                coord3 = palp(
                    self.data_dict[symbol_str(expr.c1)],
                    angle,
                    self.get_len(symbol_str(expr.c2), target),
                    self.data_dict[symbol_str(expr.c2)],
                    expr.op
                )
            else:
                return False
            if coord3.is_nan():
                return False
            t = self.mapping_r[target]
            vpoint = self.vpoints[t]
            self.data_dict[target] = coord3
            if vpoint.type == VJoint.R:
                self.data_dict[t, -1] = coord3
            else:
                self.data_dict[t, -1] = (vpoint.c[0, 0], vpoint.c[0, 1])
                self.data_dict[t, -2] = coord3
        if not self.bfgs_mode:
            return True
        # Add coordinate of known points
        p_data_dict = {}
        for i in range(len(self.vpoints)):
            # {1: 'A'} vs {'A': (10., 20.)}
            if self.mapping[i] in self.data_dict:
                p_data_dict[i] = self.data_dict[self.mapping[i]]
        # Add specified link lengths
        for k, v in self.mapping.items():
            if type(k) is frozenset:
                p_data_dict[k] = v
        # Solve
        try:
            solved_bfgs = SolverSystem(self.vpoints, {}, p_data_dict).solve()
        except ValueError:
            return False
        # Format:
        # R joint: [[p0]: (p0_x, p0_y), [p1]: (p1_x, p1_y)]
        # P or RP joint: [[p2]: ((p2_x0, p2_y0), (p2_x1, p2_y1))]
        for i in range(len(self.vpoints)):
            if self.mapping[i] in self.data_dict:
                continue
            vpoint = self.vpoints[i]
            # These points solved by Sketch Solve
            if vpoint.type == VJoint.R:
                self.data_dict[i, -1] = Coord.__new__(
                    Coord,
                    solved_bfgs[i][0],
                    solved_bfgs[i][1]
                )
            else:
                self.data_dict[i, -1] = Coord.__new__(
                    Coord,
                    solved_bfgs[i][0][0],
                    solved_bfgs[i][0][1]
                )
                self.data_dict[i, -2] = Coord.__new__(
                    Coord,
                    solved_bfgs[i][1][0],
                    solved_bfgs[i][1][1]
                )
        return True

    cdef double fitness(self, double[:] v):
        """Return fitness from chromosome.

        + Coordinates of fixed pivots. [0:self.l_base]
            [(xn, yn), ...]
        + Length and the angles of the links. [self.l_base]
        + Angle respect to the target points.
        """
        cdef int index = 0
        cdef int a_index = 0
        for m in self.mapping_list:
            if type(m) is int:
                (<VPoint>self.vpoints[m]).locate(v[index], v[index + 1])
                index += 2
            elif type(m) is str and m.startswith('A'):
                self.polar_angles[a_index] = v[index]
                a_index += 1
                index += 1
            else:
                self.mapping[m] = v[index]
                index += 1
        cdef double[:, :] angles = zeros((self.input_count, self.target_count),
                                         dtype=np_float)
        for index in range(self.input_count):
            a_index = index + self.l_base
            angles[index, :] = v[a_index:a_index + self.target_count]
        cdef double fitness = 0
        cdef int node
        target = {n: [] for n in self.target}
        for index in range(self.target_count):
            if not self.solve(angles[:, index]):
                # Punishment
                return HUGE_VAL
            for node in self.target:
                target[node].append(self.data_dict[node, -1])
        cdef Coord[:] path1, path2
        for node in self.target:
            path1 = array(target[node], dtype=object)
            if self.shape_only:
                _norm(path1, 1)
            path2 = self.target[node]
            for index in range(self.target_count):
                fitness += (<Coord>path1[index]).distance(path2[index])
        return fitness

    cpdef object result(self, double[:] v):
        """Input a generic data (variable array), return the mechanism
        expression.
        """
        cdef int index = 0
        cdef int a_index = 0
        cdef VPoint vpoint
        for m in self.mapping_list:
            if type(m) is int:
                vpoint = self.vpoints[m]
                vpoint.locate(v[index], v[index + 1])
                index += 2
            elif type(m) is str and m.startswith('A'):
                self.polar_angles[a_index] = v[index]
                a_index += 1
                index += 1
            else:
                self.mapping[m] = v[index]
                index += 1
        self.solve(v[self.l_base:self.l_base + self.target_count])
        expressions = []
        cdef int i
        cdef double x1, y1, x2, y2
        cdef Coord coord1, coord2
        for i in range(len(self.vpoints)):
            vpoint = self.vpoints[i]
            coord1 = self.data_dict[i, -1]
            vpoint.locate(coord1.x, coord1.y)
            if vpoint.type != VJoint.R:
                coord2 = self.data_dict[i, -2]
                vpoint.move((coord1.x, coord1.y), (coord2.x, coord2.y))
            expressions.append(vpoint.expr())
        return "M[" + ", ".join(expressions) + "]"
