# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, FrozenSet, Mapping, Optional, Union
from .expression import VPoint, Coord

_Coord = Tuple[float, float]
_Inputs = Mapping[Tuple[int, int], float]
_PointPair = Union[int, Tuple[int, int]]

class SolverSystem:
    def __init__(
        self,
        vpoints: Sequence[VPoint],
        inputs: Optional[Mapping[Tuple[int, int], float]] = None,
        data_dict: Optional[Mapping[_PointPair, Union[Coord, float]]] = None
    ):
        """The expression `vpoints` solver function of BFGS method by giving
        the input pairs `inputs` and link length `data_dict` requirements.

        !!! note
            The format of input pairs:

            + Revolut joints: `{(base, driver): angle}`
            + Slider joints: `{(base, base): offset}`

        The format of `data_dict`:

        + Specific coordinates: Dict\[int, List\[Coord]]
        + Specific link length: Dict\[Tuple\[int, int], float]

        The `data_dict` parameter will reformat its keys into `frozenset` type.
        """
        ...

    def same_points(self, vpoints_: Sequence[VPoint]) -> bool:
        ...

    def show_inputs(self) -> FrozenSet[Tuple[int, int]]:
        ...

    def show_data(self) -> FrozenSet[_PointPair]:
        ...

    def set_inputs(self, inputs: Mapping[Tuple[int, int], float]) -> None:
        ...

    def set_data(self, data_dict: Union[_Inputs, Mapping[int, Coord]]) -> None:
        ...

    def solve(self) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
        ...
