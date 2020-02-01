# -*- coding: utf-8 -*-

from typing import Tuple, List, Iterable, Sequence, Callable, Optional, ClassVar
from enum import IntEnum, auto

# Color type must be a RGB data
_Color = Tuple[int, int, int]
_Coord = Tuple[float, float]

def get_vlinks(vpoints: Iterable[VPoint]) -> List[VLink]:
    """Get VLinks from a list of VPoint `vpoints`."""
    ...

class Coordinate:

    """A data class used to store coordinates."""

    x: float
    y: float

    def __init__(self, x: float, y: float) -> None:
        """The constructor of Coordinate class."""
        ...

    def distance(self, p: Coordinate) -> float:
        """Return the distance between two coordinates."""
        ...

    def is_nan(self) -> bool:
        """Return True if the coordinate value is not a number."""
        ...

    def __repr__(self) -> str:
        ...

class VJoint(IntEnum):
    """Enumeration values of Joint types."""
    R = auto()
    P = auto()
    RP = auto()

class VPoint:

    """Mechanism expression class."""

    links: Sequence[str]
    c: Tuple[_Coord, _Coord]
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
        """A fast constructor of revolute joints."""
        ...

    @staticmethod
    def slider_joint(links: Iterable[str], type_int: VJoint, angle: float, x: float, y: float) -> VPoint:
        """A fast constructor of slider joints."""
        ...

    def copy(self) -> VPoint:
        """The copy method of the VPoint object."""
        ...

    @property
    def cx(self) -> float:
        """X value of current coordinate.
        If it's slider, the pin coordinate will be returned.
        """
        return 0.

    @property
    def cy(self) -> float:
        """Y value of current coordinate.
        If it's slider, the pin coordinate will be returned.
        """
        return 0.

    def set_links(self, links: Iterable[str]) -> None:
        """The update function of links attribute."""
        ...

    def replace_link(self, link1: str, link2: str) -> None:
        """Replace the value in links attribute."""
        ...

    def move(self, c1: Tuple[float, float], c2: Optional[Tuple[float, float]] = None) -> None:
        """The update function of current coordinate(s).
        The 2nd placement is the pin coordinate of slider joints.

        If there is only one argument for a slider joint,
        the slot and pin coordinates will be set to the same position.
        """
        ...

    def locate(self, x: float, y: float) -> None:
        """The update function of original coordinate.
        It will call `self.move((x, y))` after set the position.
        """
        ...

    def rotate(self, angle: float) -> None:
        """The update function of angle attribute."""
        ...

    def set_offset(self, offset: float) -> None:
        """The update function of slider offset.
        It will also enable offset value after called.
        """
        ...

    def disable_offset(self) -> None:
        """Disable offset setting of the joint."""
        ...

    def distance(self, p: VPoint) -> float:
        """Return the distance between two VPoint objects."""
        ...

    def has_offset(self) -> bool:
        """Return True if the offset setting is enabled."""
        ...

    def offset(self) -> float:
        """Return the offset constraint value of the joint."""
        ...

    def true_offset(self) -> float:
        """Return the current offset value of the joint."""
        ...

    def slope_angle(self, p: VPoint, num1: int = 2, num2: int = 2) -> float:
        """Return the value `hypot(p_x - m_x, p_y - m_y)`,
        where `m_x`, `m_y` is the value of the joint,
        and `p_x`, `p_y` is the value of `p`.

        The option `num1` and `num2` is the position of current coordinate attribute.
        """
        ...

    def grounded(self) -> bool:
        """Return True if the joint pin is connected to ground link."""
        ...

    def pin_grounded(self) -> bool:
        """Return True if the point is at the same link."""
        ...

    def same_link(self, p: VPoint) -> bool:
        """Return True if the point is at the same link."""
        ...

    def no_link(self) -> bool:
        """Return True if there is no any link in links attribute."""
        ...

    def is_slot_link(self, link_name: str) -> bool:
        """Return True if the slot is on the link `link_name`."""
        ...

    def expr(self) -> str:
        """Return the literal mechanism expression of the joint."""
        ...

    def __getitem__(self, i: int) -> float:
        ...

    def __repr__(self) -> str:
        ...

class VLink:

    """Mechanism expression class in link's view."""

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
        """The update function of points attribute."""
        ...

    def __contains__(self, point: int) -> bool:
        ...

    def __repr__(self) -> str:
        ...
