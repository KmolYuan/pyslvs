# -*- coding: utf-8 -*-

"""All collections of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import Tuple, Sequence, Dict, TypedDict, Iterator, Optional
from copy import deepcopy


class Collection(TypedDict, total=False):
    expression: str
    input: Sequence[Tuple[Tuple[int, int], Sequence[float]]]
    graph: Sequence[Tuple[int, int]]
    placement: Dict[int, Optional[Tuple[float, float, float]]]
    target: Dict[int, Optional[Sequence[Tuple[float, float]]]]
    cus: Dict[int, int]
    same: Dict[int, int]


_collection_list: Dict[str, Collection] = {
    "Four bar linkage mechanism": {
        'expression':
            "M["
            "J[R, P[0.0, 0.0], L[ground, L1]],"
            "J[R, P[90.0, 0.0], L[ground, L2]],"
            "J[R, P[12.92, 32.53], L[L1, L3]],"
            "J[R, P[73.28, 67.97], L[L2, L3]],"
            "J[R, P[33.3, 66.95], L[L3]],"
            "]",
        'input': [((0, 2), [0, 360])],
        'graph': ((0, 1), (0, 2), (1, 3), (2, 3)),
        'placement': {0: None, 1: None},
        'target': {4: None},
        'cus': {4: 3},
        'same': {},
    },
    "Six bar linkage mechanism": {
        'expression':
            "M["
            "J[R, P[5.8511, -103.2831], L[ground, L2]],"
            "J[R, P[-71.3292, -109.7644], L[ground, L3]],"
            "J[R, P[77.8903, -110.263], L[ground, L5]],"
            "J[R, P[23.7994, 20.3606], L[L1, L2]],"
            "J[R, P[-33.1974, 90.658], L[L1, L4]],"
            "J[R, P[112.2951, 18.3663], L[L1, L5]],"
            "J[R, P[-117.0519, -18.5671], L[L3, L4]]"
            "]",
        'graph': ((0, 2), (2, 1), (0, 3), (3, 4), (4, 1), (0, 5), (5, 1)),
        'placement': {0: None, 1: None, 2: None},
        'target': {4: None},
        'cus': {},
        'input': [((1, 6), [0, 360])],
        'same': {},
    },
    "Eight bar linkage mechanism": {
        'expression':
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
        'input': [((0, 3), [0, 360])],
        'graph': ((0, 1), (0, 4), (0, 5), (1, 2), (1, 3), (2, 4), (3, 5),
                  (3, 7), (4, 6), (6, 7)),
        'placement': {0: None, 1: None},
        'target': {10: None},
        'cus': {10: 7},
        'same': {2: 1, 4: 3, 7: 6},
    },
    "Ball lifter linkage mechanism": {
        'expression':
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
        'input': [((0, 5), [0, 360])],
        'graph': ((0, 1), (0, 4), (0, 6), (0, 7), (0, 9), (1, 2), (1, 3),
                  (2, 4), (2, 5), (3, 7), (3, 8), (5, 6), (8, 9)),
        'placement': {0: None, 1: None, 2: None, 3: None, 4: None},
        'target': {13: None, 14: None},
        'cus': {13: 5, 14: 8},
        'same': {6: 5},
    },
}


def collection_list(key: str) -> Collection:
    """The example data of collections.

    The format of each configuration is:

    + `expression`: Mechanism expression of the structure.
        + type: str
    + `input`: Input pairs.
        + type: Sequence[Tuple[int, int]]
    + `graph`: The generalized chain graph in edge set.
        + type: Sequence[Tuple[int, int]]
    + `placement`: The grounded joints setting. (`x`, `y`, `r`)
        + type: Dict[int, Optional[Tuple[float, float, float]]]
    + `target`: The target joints settings.
        + type: Dict[int, Optional[Sequence[Tuple[float, float]]]]
    + `cus`: The custom joints on specific link. (link number correspond to
        the graph expression.)
        + type: Dict[int, int]
    + `same`: The multiple joints setting.
        + type: Dict[int, int]
    """
    return deepcopy(_collection_list[key])


def all_collections() -> Iterator[str]:
    """Get all collection names."""
    yield from sorted(_collection_list)
