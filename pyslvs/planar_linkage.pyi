# -*- coding: utf-8 -*-

from typing import Tuple, List, Dict, Iterable, Any
from numpy import ndarray, double
from .metaheuristics import ObjFunc

_Coord = Tuple[float, float]

def norm_path(path: Iterable[_Coord], scale: float = 1) -> List[_Coord]:
    ...

def curvature(path: Iterable[_Coord]) -> ndarray:
    ...

def derivative(path: ndarray) -> ndarray:
    ...

def path_signature(k: ndarray, maximum: float = 100) -> ndarray:
    ...

def cross_correlation(p1: ndarray, p2: ndarray, t: float = 0.1) -> ndarray:
    ...

class FMatch(ObjFunc[str]):

    def __init__(self, mech: Dict[str, Any]):
        """The constructor of objective object.

        Options of `mech_params`:

        + `Expression`: The mechanism expression of the structure.
            + type: List\[[VPoint]]
        + `input`: Input pairs.
            + type: List[Tuple[int, int]]
        + `Placement`: The grounded joints setting. (`x`, `y`, `r`)
            + type: Dict[int, Tuple[float, float, float]]
        + `Target`: The target path.
            + type: Dict[int, Sequence[Tuple[float, float]]]
        + `same`: Multiple joint setting. The joints are according to [`edges_view`](#edges_view).
            + type: Dict[int, int]
        + `upper`: The upper setting of variables, the length must same as variable array.
            + type: List[float]
        + `lower`: The lower setting of variables, the length must same as variable array.
            + type: List[float]
        + `shape_only`: Compare paths by shape only.
            + type: bool

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
