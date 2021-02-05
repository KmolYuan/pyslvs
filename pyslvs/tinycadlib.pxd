# -*- coding: utf-8 -*-
# cython: language_level=3

"""Sharing position analysis function.

author: Yuan Chang
copyright: Copyright (C) 2016-2021
license: AGPL
email: pyslvs@gmail.com
"""

from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.vector cimport vector
from .topo_config cimport EStack

cdef extern from "tinycadlib/solver.h" nogil:
    struct SwappablePair:
        int first, second

    enum Label:
        P_LABEL
        L_LABEL
        I_LABEL
        A_LABEL
        S_LABEL

    enum Func:
        PXY
        PPP
        PLA
        PLAP
        PLLP
        PLPP
        PALP

    struct CCoord:
        double x, y

    CCoord cpxy(CCoord c1, double x, double y)
    CCoord cppp(CCoord c1, CCoord c2, CCoord c3)
    CCoord cplap(CCoord c1, double d0, double a0, CCoord c2, bint inverse)
    CCoord cpllp(CCoord c1, double d0, double d1, CCoord c2, bint inverse)
    CCoord cplpp(CCoord c1, double d0, CCoord c2, CCoord c3, bint inverse)
    CCoord cpalp(CCoord c1, double a0, double d0, CCoord c2, bint inverse)

    ctypedef pair[int, int] Sym

    struct Expr:
        bint op
        Func func
        Sym v1, v2, c1, c2, c3, target

    cppclass ExprSolver:
        map[Sym, CCoord] joint_pos

        ExprSolver()
        ExprSolver(vector[Expr] stack, map[Sym, CCoord] joint_pos,
                   map[Sym, double] param)
        bint solve()

cdef bint preprocessing(EStack exprs, object vpoints, object inputs,
                        map[Sym, CCoord] & joint_pos,
                        map[SwappablePair, double] & link_len,
                        map[Sym, double] & param)
cpdef list expr_solving(EStack exprs, object vpoints, object inputs=*)
cdef (bint, map[Sym, CCoord]) quick_solve(vector[Expr] stack,
                                          map[Sym, CCoord] joint_pos,
                                          map[Sym, double] param) nogil
cdef double[:, :, :] c_uniform_path(double[:, :] dimension, int n) nogil
cpdef object uniform_expr(double[:] v)
