# -*- coding: utf-8 -*-

from typing import Tuple, Iterable
from numpy import ndarray, double
from pyslvs.metaheuristics import ObjFunc
from .utility import FConfig

_Coord = Tuple[float, float]

def norm_path(path: Iterable[_Coord], scale: float = 1) -> ndarray:
    ...

def curvature(path: Iterable[_Coord]) -> ndarray:
    ...

def derivative(path: ndarray) -> ndarray:
    ...

def path_signature(k: ndarray, maximum: float = 100) -> ndarray:
    ...

def cross_correlation(p1: ndarray, p2: ndarray, t: float = 0.1) -> ndarray:
    ...


class FPlanar(ObjFunc[str]):
    callback: int

    def __init__(self, mech: FConfig):
        """The constructor of objective object.

        Options:
        + `expression`: The mechanism expression of the structure.
        + `input`: Input pairs.
        + `placement`: The grounded joints setting. (`x`, `y`, `r`)
        + `target`: The target path.
        + `same`: Multiple joint setting. The joints are according to [`edges_view`](#edges_view).
        + `upper`: The upper setting of variables, the length must same as variable array.
        + `lower`: The lower setting of variables, the length must same as variable array.
        + `shape_only`: Compare paths by shape only.

        Variable array:

        | | Placement | Link length | Inputs |
        |:---:|:-----:|:-----------:|:------:|
        | `v =` | `x0`, `y0`, ... | `l0`, `l1`, ... | `a00`, `a01`, ..., `a10`, `a11`, ... |

        In 1D array: `v = [x0, y0, ..., l0, l1, ..., a00, a01, ..., a10, a11, ...]`
        """
        ...

    def fitness(self, v: ndarray) -> double:
        """The fitness is the error between target path and self.

        Chromosome format: (decided by upper and lower)

        v: `[Ax, Ay, Dx, Dy, ..., L0, L1, ..., A00, A01, ..., A10, A11, ...]`
        """
        ...

    def is_two_kernel(self) -> bool:
        ...

    def result(self, v: ndarray) -> str:
        ...
