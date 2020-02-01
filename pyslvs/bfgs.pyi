# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, FrozenSet, Dict, Optional, Union
from .expression import VPoint, Coordinate

_Coord = Tuple[float, float]
_Inputs = Dict[Tuple[int, int], float]

class SolverSystem:

    """Sketch Solve solver.

    !!! note
        The object attributes of such type are unable to access.
    """

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
        """Return True if two expressions are same."""
        ...

    def show_inputs(self) -> FrozenSet[Tuple[int, int]]:
        """Show the current input pairs keys from original constructor."""
        ...

    def show_data(self) -> FrozenSet[Union[int, Tuple[int, int]]]:
        """Show the current keys of `data_dict` parameter from original constructor."""
        ...

    def set_inputs(self, inputs: Dict[Tuple[int, int], float]) -> None:
        """Set the values of `inputs` parameter from original constructor.
        Two groups of `dict` keys must be the same or subset.
        """
        ...

    def set_data(self, data_dict: Union[_Inputs, Dict[int, Coordinate]]) -> None:
        """Set the values of `data_dict` parameter from original constructor.
        Two groups of `dict` keys must be the same or subset.
        """
        ...

    def solve(self) -> List[Union[_Coord, Tuple[_Coord, _Coord]]]:
        """Solve the conditions and return the result, raise ValueError if not succeeded.
        The joint position will returned by its index correspondingly.

        + Revolut joints: Tuple[float, float]
        + Slider joints: Tuple[Tuple[float, float], Tuple[float, float]]
        """
        ...
