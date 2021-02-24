# -*- coding: utf-8 -*-

from typing import Iterable, Tuple
from numpy import double, ndarray
from pyslvs.metaheuristics import ObjFunc
from .utility import NConfig

_Coord = Tuple[float, float]

def norm_pca(path: Iterable[_Coord]) -> ndarray:
    ...


class NPlanar(ObjFunc[str]):
    def __init__(self, mech: NConfig):
        """The constructor of objective object.

        Options:
        + `target`: The target path.
        """

    def fitness(self, v: ndarray) -> double:
        ...

    def result(self, v: ndarray) -> str:
        ...
