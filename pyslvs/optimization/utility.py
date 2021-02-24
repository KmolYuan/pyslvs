# -*- coding: utf-8 -*-

"""Optimization utilities."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import TypedDict, Sequence, Tuple, Dict
from pyslvs.expression import VPoint

_Path = Sequence[Tuple[float, float]]


class FConfig(TypedDict, total=False):
    expression: Sequence[VPoint]
    input: Sequence[Tuple[Tuple[int, int], Sequence[float]]]
    placement: Dict[int, Tuple[float, float, float]]
    target: Dict[int, _Path]
    same: Dict[int, int]
    upper: float
    lower: float
    shape_only: bool


class NConfig(TypedDict, total=False):
    target: _Path
