# -*- coding: utf-8 -*-

from typing import Tuple, List, Sequence, Dict, Union, Optional
from .triangulation import ExpressionStack
from .expression import VPoint, Coordinate

TuplePoint = Tuple[float, float]

def plap(
    c1: Coordinate,
    d0: float,
    a0: float,
    c2: Optional[Coordinate] = None,
    inverse: bool = False
) -> Coordinate:
    """The PLAP function requires two points, one distance and one angle,
    obtained the position of third point. The unit of `a0` is degree.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `d0` correspond to "L0", `a0` correspond to "beta", `return` correspond to "C".
    If `c2` is not given, "alpha" will be set to zero.

    ![PLAP](img/PLAP.png)

    Set `inverse` option to `True` can make `a0` value as negative.
    """
    ...

def pllp(
    c1: Coordinate,
    d0: float,
    d1: float,
    c2: Coordinate,
    inverse: bool = False
) -> Coordinate:
    """The PLLP function requires two points and two distances, obtained the position of third point.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `d0` correspond to "L0", `d1` correspond to "L1", `return` correspond to "C".

    ![PLLP](img/PLLP.png)

    Set `inverse` option to `True` can make the result upside down.
    """
    ...

def plpp(
    c1: Coordinate,
    d0: float,
    c2: Coordinate,
    c3: Coordinate,
    inverse: bool = False
) -> Coordinate:
    """The PLLP function requires three points and one distance, obtained the position of fourth point.

    In the following picture, `c1` correspond to "A", `c2` correspond to "B",
    `c3` correspond to "C", `d0` correspond to "L0", `return` correspond to "D".

    ![PLPP](img/PLPP.png)

    Set `inverse` option to `True` can make the result to the another side
    between `c1` and line `c2` `c3`.
    """
    ...

def pxy(c1: Coordinate, x: float, y: float) -> Coordinate:
    """The PXY function requires one point and offset values, obtained the position of second point.

    In the following picture, `c1` correspond to "A", `d0` correspond to "X",
    `d1` correspond to "Y", `return` correspond to "B", the sign of value are
    correspond to coordinate system.

    ![PXY](img/PXY.png)
    """
    ...

def vpoint_dof(vpoints: Sequence[VPoint]) -> int:
    """Return the DOF of the mechanism expression `vpoints`."""
    ...

def expr_parser(exprs: Sequence[Tuple[str, ...]], data_dict: Dict[str, float]) -> None:
    """Solve and update information of the triangle expression `exprs` to `data_dict`.
    The argument `exprs` can be obtained by [`vpoints_configure`](#vpoints_configure)
    and [`ExpressionStack.as_list()`](#expressionstackas_list) method.

    This function is already included in [`expr_solving`](#expr_solving),
    not recommended for direct use.
    """
    ...

def data_collecting(
    exprs: ExpressionStack,
    mapping: Dict[int, str],
    vpoints_: Sequence[VPoint],
) -> Tuple[Dict[str, Union[Coordinate, float]], int]:
    """Data transform function of Triangular method.
    The triangle expression stack `expr` is generated from [`vpoints_configure`](#vpoints_configure).
    The information data `mapping` map the symbols to the indicator of `vpoints_`.

    This function is already included in [`expr_solving`](#expr_solving),
    not recommended for direct use.
    """
    ...

def expr_solving(
    exprs: ExpressionStack,
    mapping: Dict[Union[int, Tuple[int, int]], Union[str, float]],
    vpoints: Sequence[VPoint],
    angles: Optional[Sequence[float]] = None
) -> List[Union[TuplePoint, Tuple[TuplePoint, TuplePoint]]]:
    """Solver function of Triangular method and BFGS method, for mechanism expression `vpoints`.

    The triangle expression stack `expr` is generated from [`vpoints_configure`](#vpoints_configure).

    The information data `mapping` map the symbols to the indicator of `vpoints`,
    additionally has a same format as argument `data_dict` in [SolverSystem].

    Solver function will not handle slider input pairs in argument `angles`,
    which is only support revolute joints. In another way, the slider input pairs
    can be set by [`VPoint.disable_offset()`](#vpointdisable_offset) method.
    """
    ...
