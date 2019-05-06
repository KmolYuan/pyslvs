# -*- coding: utf-8 -*-

"""All collections of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"


collection_list = {
    "Four bar linkage mechanism": {
        'Expression':
            "M["
            "J[R, P[-70, -70], L[ground, L1]],"
            "J[R, P[70, -70], L[ground, L2]],"
            "J[R, P[-70, 12.5], L[L1, L3]],"
            "J[R, P[70, 12.5], L[L2, L3]],"
            "J[R, P[0, 63.5], L[L3]]"
            "]",
        'input': [(0, 2)],
        'Graph': ((0, 1), (0, 2), (1, 3), (2, 3)),
        'Placement': {0: None, 1: None},
        'Target': {4: None},
        'cus': {4: 3},
        'same': {},
    },

    "Six bar linkage mechanism": {
        'Expression':
            "M["
            "J[R, P[5.8511, -103.2831], L[ground, L2]], "
            "J[R, P[-71.3292, -109.7644], L[ground, L3]], "
            "J[R, P[77.8903, -110.263], L[ground, L5]], "
            "J[R, P[23.7994, 20.3606], L[L1, L2]], "
            "J[R, P[-33.1974, 90.658], L[L1, L4]], "
            "J[R, P[112.2951, 18.3663], L[L1, L5]], "
            "J[R, P[-117.0519, -18.5671], L[L3, L4]]"
            "]",
        'Graph': ((0, 2), (2, 1), (0, 3), (3, 4), (4, 1), (0, 5), (5, 1)),
        'Placement': {0: None, 1: None, 2: None},
        'Target': {4: None},
        'cus': {},
        'input': [(1, 6)],
        'same': {},
    },

    "Eight bar linkage mechanism": {
        'Expression':
            "M["
            "J[R, P[30.5, 10.5], L[ground, L1]],"
            "J[R, P[-14.5, 10.5], L[ground, L4, L5]],"
            "J[R, P[81.5, 60.5], L[L1, L2, L3]],"
            "J[R, P[-31.5, 86.5], L[L2, L4]],"
            "J[R, P[41.5, -38.5], L[L3, L5, L7]],"
            "J[R, P[-85.5, 9.5], L[L4, L6]],"
            "J[R, P[-37.5, -48.5], L[L6, L7]],"
            "J[R, P[35.5, -107.5], L[L7]]"
            "]",
        'input': [(0, 3)],
        'Graph': (
            (0, 1),
            (0, 4),
            (0, 5),
            (1, 2),
            (1, 3),
            (2, 4),
            (3, 5),
            (3, 7),
            (4, 6),
            (6, 7),
        ),
        'Placement': {0: None, 1: None},
        'Target': {10: None},
        'cus': {10: 7},
        'same': {2: 1, 4: 3, 7: 6},
    },

    "Ball lifter linkage mechanism": {
        'Expression':
            "M["
            "J[R, P[36.5, -59.5], L[ground, L1]],"
            "J[R, P[10, -94.12], L[ground, L4]],"
            "J[R, P[-28.5, -93.5], L[ground, L6]],"
            "J[R, P[102.5, -43.5], L[ground, L7]],"
            "J[R, P[77.5, -74.5], L[ground, L9]],"
            "J[R, P[28.82, -22.35], L[L1, L2, L3]],"
            "J[R, P[-18.5, -44.5], L[L2, L4]],"
            "J[R, P[-75.5, -59.5], L[L2, L5]],"
            "J[R, P[56.5, 29.5], L[L3, L7]],"
            "J[R, P[68.5, 71.5], L[L8, L3]],"
            "J[R, P[-47.06, -28.24], L[L5, L6]],"
            "J[R, P[107.5, 42.5], L[L8, L9]],"
            "J[R, P[-109.41, -49.41], L[L5]],"
            "J[R, P[44.12, 107.65], L[L8]]"
            "]",
        'input': [(0, 5)],
        'Graph': (
            (0, 1),
            (0, 4),
            (0, 6),
            (0, 7),
            (0, 9),
            (1, 2),
            (1, 3),
            (2, 4),
            (2, 5),
            (3, 7),
            (3, 8),
            (5, 6),
            (8, 9),
        ),
        'Placement': {0: None, 1: None, 2: None, 3: None, 4: None},
        'Target': {13: None, 14: None},
        'cus': {13: 5, 14: 8},
        'same': {6: 5},
    },
}
