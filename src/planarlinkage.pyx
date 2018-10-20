# -*- coding: utf-8 -*-
# cython: language_level=3

"""The callable classes of the validation in algorithm.

author: Yuan Chang
copyright: Copyright (C) 2016-2018
license: AGPL
email: pyslvs@gmail.com
"""

cimport cython
from numpy import (
    array as np_array,
    object as np_object,
    float32 as np_float32,
)
from numpy cimport ndarray
from verify cimport Verification
from tinycadlib cimport (
    Coordinate,
    radians,
    PLAP,
    PLLP,
    PLPP,
    PXY,
    legal_crank,
    str_between,
    str_before,
)

# A large fitness. Infinity can not used in chart.
cdef double FAILURE = 9487945


@cython.final
cdef class Planar(Verification):

    """This class is used to verified kinematics of the linkage mechanism."""

    cdef int target_count, var_count
    cdef list link_list, driver_list, follower_list
    cdef ndarray constraints, target_names, exprs, target, upper, lower

    def __cinit__(self, mech_params: dict):
        """mech_params = {
            'Driver': {'pt': (x, y, r)},
            'Follower': {'pt': (x, y, r)},
            'Target': {'pt': [(x0, y0), (x1, y1), ...]},
            'constraints': [('pt', 'pt', 'pt', 'pt')],
            'Expression': str,
            'upper': ndarray[double, ndim=1],
            'lower': ndarray[double, ndim=1],
        }
        """
        cdef set check_set = set(map(len, mech_params['Target'].values()))
        if len(check_set) != 1:
            raise ValueError("Target path should be in the same size.")
        # Counting how many action to satisfied require point.
        self.target_count = check_set.pop()
        # Target points
        # [Coordinate(x0, y0), Coordinate(x1, y1), Coordinate(x2, y2), ...]
        cdef int i = 0
        cdef int target_count = len(mech_params['Target'])
        self.target_names = ndarray(target_count, dtype=np_object)
        self.target = ndarray(target_count, dtype=np_object)
        cdef double x, y
        cdef str name
        cdef list target
        for name, target in mech_params['Target'].items():
            self.target_names[i] = name
            self.target[i] = tuple(Coordinate(x, y) for x, y in target)
            i += 1

        # Constraint
        self.constraints = np_array(mech_params['constraints'])

        # Expression
        # ['A', 'B', 'C', 'D', 'E', 'L0', 'L1', 'L2', 'L3', 'L4', 'a0']
        cdef tuple exprs = tuple(mech_params['Expression'].split(';'))

        """
        link_list: L0, L1, L2, L3, ...
        driver_list: The name of the point in "self.Driver" (Sorted).
        follower_list: The name of the point in "self.Follower" (Sorted).
        exprs:
            {'relate': 'PLAP', 'target': 'B', 'params': ['A', 'L0', 'a0', 'D']},
            {'relate': 'PLLP', 'target': 'C', 'params': ['B', 'L1', 'L2', 'D']}, ...
        Expression: PLAP[A,L0,a0,D](B);PLLP[B,L1,L2,D](C);PLLP[B,L3,L4,C](E)
        """
        self.link_list = []
        self.driver_list = []
        self.follower_list = []
        self.exprs = ndarray(len(exprs), dtype=np_object)
        cdef str expr, params, p
        for i, expr in enumerate(exprs):
            params = str_between(expr, '[', ']')
            self.exprs[i] = (
                # [0]: relate
                str_before(expr, '['),
                # [1]: target
                str_between(expr, '(', ')'),
                # [2]: params
                params.split(','),
            )
            # Parameters should exist in solver expression.
            for p in params.split(','):
                if p.startswith('L') and len(p) > 1:
                    self.link_list.append(p)
                if (p in mech_params['Driver']) and (p not in self.driver_list):
                    self.driver_list.append(p)
                if (p in mech_params['Follower']) and (p not in self.follower_list):
                    self.follower_list.append(p)

        """Limitations
        
        self.var_count = length before matching angles.
        
        The data format will same as chromosome.
        [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        Input format of upper and lower:
        [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A0, A1, ...]
        + The input angles are not duplicated.
        + The number of angles will same as length of driver list.
        """
        cdef int link_count = len(self.link_list)
        # The number of all variables (except angles).
        self.var_count = 2 * len(self.driver_list) + 2 * len(self.follower_list) + link_count

        cdef dict tmp_dict = {}
        tmp_dict.update(mech_params['Driver'])
        tmp_dict.update(mech_params['Follower'])

        cdef list tmp_list = []

        # upper
        for name in self.driver_list + self.follower_list:
            for i in range(2):
                tmp_list.append(tmp_dict[name][i] + tmp_dict[name][2]/2)
        tmp_list.extend(mech_params['upper'][:link_count])
        for i in range(len(self.driver_list)):
            tmp_list.extend([mech_params['upper'][link_count + i]] * self.target_count)
        self.upper = np_array(tmp_list, dtype=np_float32)

        tmp_list.clear()

        # lower
        for name in self.driver_list + self.follower_list:
            for i in range(2):
                tmp_list.append(tmp_dict[name][i] - tmp_dict[name][2]/2)
        tmp_list.extend(mech_params['lower'][:link_count])
        for i in range(len(self.driver_list)):
            tmp_list.extend([mech_params['lower'][link_count + i]] * self.target_count)
        self.lower = np_array(tmp_list, dtype=np_float32)

        # Swap sorting.
        for i in range(len(self.upper)):
            if self.upper[i] < self.lower[i]:
                self.upper[i], self.lower[i] = self.lower[i], self.upper[i]

    cdef ndarray[double, ndim=1] get_upper(self):
        return self.upper

    cdef ndarray[double, ndim=1] get_lower(self):
        return self.lower

    cdef int get_nParm(self):
        return len(self.upper)

    cdef dict get_data_dict(self, ndarray[double, ndim=1] v):
        """Create and return data dict."""
        cdef str name
        cdef dict tmp_dict = {}
        cdef int vi = 0
        # driver and follower
        for name in self.driver_list + self.follower_list:
            tmp_dict[name] = Coordinate(v[vi], v[vi+1])
            vi += 2
        # links
        for name in self.link_list:
            tmp_dict[name] = v[vi]
            vi += 1
        return tmp_dict

    cdef inline ndarray get_path_array(self):
        """Create and return path array."""
        cdef ndarray path = ndarray(len(self.target_names), dtype=np_object)
        cdef int i
        for i in range(len(self.target_names)):
            path[i] = []
        return path

    cdef Coordinate from_formula(self, tuple expr, dict data_dict):
        """Formulas using PLAP and PLLP."""
        cdef str fun = expr[0]
        cdef str p
        cdef list params = []
        for p in expr[2]:
            if p == 'T':
                params.append(True)
            elif p == 'F':
                params.append(False)
            else:
                params.append(data_dict[p])
        cdef int params_count = len(params)

        # We should unpack as C++'s way.
        cdef double x = float('nan')
        cdef double y = x
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
        return Coordinate(x, y)

    cdef double run(self, ndarray[double, ndim=1] v) except *:
        """Chromosome format: (decided by upper and lower)
        
        v: [Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]
        target: a list of target. [(1,5), (2,5), (3,5)]
        target_count: length of target.
        vars: mechanism variables count.
        """
        # all variable
        cdef dict test_dict = self.get_data_dict(v)
        cdef ndarray path = self.get_path_array()
        # calculate the target point, and sum all error.
        # My fitness
        cdef double fitness = 0
        cdef int i, j, k
        cdef str name
        cdef tuple e
        cdef Coordinate target_coord
        for i in range(self.target_count):
            # a0: random angle to generate target point.
            # match to path points.
            for j in range(len(self.driver_list)):
                test_dict[f'a{j}'] = radians(v[self.var_count + i * len(self.driver_list) + j])
            for e in self.exprs:
                # target
                target_coord = self.from_formula(e, test_dict)
                if target_coord.is_nan():
                    return FAILURE
                else:
                    test_dict[e[1]] = target_coord
            for i, name in enumerate(self.target_names):
                path[i].append(test_dict[name])
        # constraint
        cdef ndarray constraint
        for constraint in self.constraints:
            if not legal_crank(
                test_dict[constraint[0]],
                test_dict[constraint[1]],
                test_dict[constraint[2]],
                test_dict[constraint[3]]
            ):
                return FAILURE
        # swap
        cdef list errors
        cdef Coordinate c
        for k in range(len(self.target_names)):
            for i in range(self.target_count):
                fitness += path[k][i].distance(self.target[k][i])
        return fitness

    cpdef object result(self, ndarray v):
        """Return the last answer."""
        cdef str k
        cdef tuple e
        cdef object value
        cdef dict final_dict = self.get_data_dict(v)
        for j in range(len(self.driver_list)):
            final_dict[f'a{j}'] = radians(v[self.var_count + j])
        for e in self.exprs:
            # target
            final_dict[e[1]] = self.from_formula(e, final_dict)
        for k, value in final_dict.items():
            if type(value) == Coordinate:
                final_dict[k] = (value.x, value.y)
        cdef list tmp_list
        for j in range(len(self.driver_list)):
            tmp_list = []
            for i in range(self.target_count):
                tmp_list.append(v[self.var_count + i * len(self.driver_list) + j])
            final_dict[f'a{j}'] = tuple(tmp_list)
        return final_dict

    def __call__(self, v: ndarray) -> double:
        """Python callable object."""
        return self.run(v)
