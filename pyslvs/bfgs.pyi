# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, FrozenSet, Dict, Optional, Union
from .expression import VPoint, Coordinate

_Coord = Tuple[float, float]
_Inputs = Dict[Tuple[int, int], float]

class SolverSystem:

    def __init__(
        self,
        vpoints: Sequence[VPoint],
        inputs: Optional[Dict[Tuple[int, int], float]] = None,
        data_dict: Optional[Dict[Union[int, Tuple[int, int]], Union[Coordinate, float]]] = None
    ) -> None:
        """The expression `vpoints` solver function of BFGS method by
        giving the input pairs `inputs` and link length `data_dict` requirements.

        !!! note
            The format of input pairs:

            + Revolut joints: `{(base, driver): angle}`
            + Slider joints: `{(base, base): offset}`

        The format of `data_dict`:

        + Specific coordinates: Dict\[int, [Coordinate]]
        + Specific link length: Dict\[Tuple[int, int], float]

        The `data_dict` parameter will reformat its keys into `frozenset` type.
        """
        ...

    def same_points(self, vpoints_: Sequence[VPoint]) -> bool:
        ...

    def show_inputs(self) -> FrozenSet[Tuple[int, int]]:
        ...

    def show_data(self) -> FrozenSet[Union[int, Tuple[int, int]]]:
        ...

    def set_inputs(self, inputs: Dict[Tuple[int, int], float]) -> None:
        ...

    def set_data(self, data_dict: Union[_Inputs, Dict[int, Coordinate]]) -> None:
        ...

    def solve(self) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
        ...
