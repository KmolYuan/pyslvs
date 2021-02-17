# -*- coding: utf-8 -*-

from typing import Iterable, Tuple
from numpy import ndarray

_Coord = Tuple[float, float]

def norm_pca(path: Iterable[_Coord]) -> ndarray:
    ...
