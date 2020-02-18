# -*- coding: utf-8 -*-

from typing import Tuple, List, Dict, Iterable, Any
from numpy import ndarray, double
from .metaheuristics import Objective


def norm_path(path: Iterable[Tuple[float, float]], scale: float = 1) -> List[Tuple[float, float]]:
    """Python wrapper of normalization function."""
    ...


class Planar(Objective[str]):

    """This class is used to verified kinematics of the linkage mechanism."""

    def __init__(self, mech: Dict[str, Any]) -> None:
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

    def get_upper(self) -> ndarray:
        """Return upper bound."""
        ...

    def get_lower(self) -> ndarray:
        """Return lower bound."""
        ...

    def fitness(self, v: ndarray) -> double:
        """The fitness is the error between target path and self."""
        ...

    def is_two_kernel(self) -> bool:
        """Input a generic data (variable array), return the mechanism expression."""
        ...

    def result(self, v: ndarray) -> str:
        """Input a generic data (variable array), return the mechanism expression."""
        ...
