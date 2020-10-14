# -*- coding: utf-8 -*-

from typing import Tuple, List, Iterable, Sequence, Callable, Optional, ClassVar
from enum import IntEnum, auto
from numpy import ndarray

# Color type must be a RGB data
_Color = Tuple[int, int, int]

def get_vlinks(vpoints: Iterable[VPoint]) -> List[VLink]:
    ...

class Coord:
    x: float
    y: float

    def __init__(self, x: float, y: float):
        """The constructor of Coordinate class."""
        ...

    def distance(self, p: Coord) -> float:
        ...

    def slope_angle(self, p: Coord) -> float:
        ...

    def is_nan(self) -> bool:
        ...

class VJoint(IntEnum):
    """Enumeration values of Joint types."""
    R = auto()
    P = auto()
    RP = auto()

class VPoint:
    links: Sequence[str]
    c: ndarray
    type: VJoint
    color: Optional[_Color]
    color_str: str
    type_str: str
    x: float
    y: float
    angle: float
    HOLDER: ClassVar[VPoint] = ...

    def __init__(
        self,
        links: Iterable[str],
        type_int: VJoint,
        angle: float,
        color_str: str,
        x: float,
        y: float,
        color_func: Optional[Callable[[str], _Color]] = None
    ):
        """The attributes will match to the object attributes of [VPoint] objects.

        Where the color function `color_func` needs to transform the color string `color_str` into RGB format.
        If color information is not needed, the `color_func` can be `None`.

        !!! note
            Some of the attributes are not work in some of the joint types.
        """
        ...

    @staticmethod
    def r_joint(links: Iterable[str], x: float, y: float) -> VPoint:
        ...

    @staticmethod
    def slider_joint(
        links: Iterable[str],
        type_int: VJoint,
        angle: float,
        x: float,
        y: float
    ) -> VPoint:
        ...

    def copy(self) -> VPoint:
        ...

    @property
    def sx(self) -> float:
        ...

    @property
    def sy(self) -> float:
        ...

    @property
    def cx(self) -> float:
        ...

    @property
    def cy(self) -> float:
        ...

    def set_links(self, links: Iterable[str]) -> None:
        ...

    def replace_link(self, link1: str, link2: str) -> None:
        ...

    def move(
        self,
        c1: Tuple[float, float],
        c2: Optional[Tuple[float, float]] = None
    ) -> None:
        ...

    def locate(self, x: float, y: float) -> None:
        ...

    def rotate(self, angle: float) -> None:
        ...

    def set_offset(self, offset: float) -> None:
        ...

    def disable_offset(self) -> None:
        ...

    def is_slider(self) -> bool:
        ...

    def distance(self, p: VPoint) -> float:
        ...

    def has_offset(self) -> bool:
        ...

    def offset(self) -> float:
        ...

    def true_offset(self) -> float:
        ...

    def slope_angle(self, p: VPoint, num1: int = 2, num2: int = 2) -> float:
        ...

    def link_pos(self, link: str) -> Coord:
        ...

    def grounded(self) -> bool:
        ...

    def pin_grounded(self) -> bool:
        ...

    def same_link(self, p: VPoint) -> bool:
        ...

    def no_link(self) -> bool:
        ...

    def is_slot_link(self, link: str) -> bool:
        ...

    def expr(self) -> str:
        ...

    def to_coord(self, ind: int) -> Coord:
        ...

    def __getitem__(self, i: int) -> float:
        ...

class VLink:
    name: str
    color_str: str
    color: Optional[_Color]
    points: Sequence[int]
    HOLDER: ClassVar[VLink] = ...
    FRAME: ClassVar[str] = ...

    def __init__(
        self,
        name: str,
        color_str: str,
        points: Iterable[int],
        color_func: Optional[Callable[[str], _Color]] = None
    ):
        """The attributes will match to the object attributes of [VLink] objects.

        Where the color function `color_func` needs to transform the color string `color_str` into RGB format.
        If color information is not needed, the `color_func` can be `None`.
        """
        ...

    def set_points(self, points: Iterable[int]) -> None:
        ...

    def points_pos(self, vpoints: Iterable[VPoint]) -> Sequence[Coord]:
        ...

    def __contains__(self, point: int) -> bool:
        ...
