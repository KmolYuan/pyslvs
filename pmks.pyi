# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    Callable,
    Optional,
)
from numpy import ndarray
from enum import IntEnum

# Color type must be a RGB data.
_Color = Tuple[int, int, int]


class VJoint(IntEnum):
    R: int
    P: int
    RP: int


class VPoint:

    """Symbol of joints."""

    links: Tuple[str, ...]
    c: ndarray
    type: int
    color: Optional[_Color]
    colorSTR: str
    typeSTR: str
    x: float
    y: float
    angle: float

    def __init__(
        self,
        links: str,
        type_int: int,
        angle: float,
        color_str: str,
        x: float,
        y: float,
        color_func: Optional[Callable[..., _Color]] = None
    ):
        ...

    @staticmethod
    def r_joint(links: str, x: float, y: float) -> VPoint:
        """Create by coordinate."""
        ...

    @staticmethod
    def slider_joint(links: str, type_int: int, angle: float, x: float, y: float) -> VPoint:
        """Create by coordinate."""
        ...

    def copy(self) -> VPoint:
        """Copy method of Python."""
        ...

    @property
    def cx(self) -> float:
        """X value of first current coordinate."""
        return 0

    @property
    def cy(self) -> float:
        """Y value of first current coordinate."""
        return 0

    def move(self, c1: Tuple[float, float], c2: Optional[Tuple[float, float]] = None):
        """Change coordinates of this point."""
        ...

    def rotate(self, angle: float):
        """Change the angle of slider slot by degrees."""
        ...

    def set_offset(self, offset: float):
        """Set slider offset."""
        ...

    def disable_offset(self):
        """Disable offset status."""
        ...

    def distance(self, p: VPoint) -> float:
        """Distance between two VPoint."""
        ...

    def has_offset(self) -> bool:
        """Return has offset."""
        ...

    def offset(self) -> float:
        """Return target offset."""
        ...

    def true_offset(self) -> float:
        """Return offset between slot and pin."""
        ...

    def slope_angle(self, p: VPoint, num1: int = 2, num2: int = 2) -> float:
        """Angle between horizontal line and two point.

        num1: me.
        num2: other side.
        [0]: base (slot) link.
        [1]: pin link.
        """
        ...

    def grounded(self) -> bool:
        """Return True if the joint is connect with the ground."""
        ...

    def pin_grounded(self) -> bool:
        """Return True if the joint has any pin connect with the ground."""
        ...

    def same_link(self, p: VPoint) -> bool:
        """Return True if the point is at the same link."""
        ...

    def no_link(self) -> bool:
        """Return True if the point has no link."""
        ...

    def is_slot_link(self, link_name: str) -> bool:
        """Return True if the link name is first link."""
        ...

    @property
    def expr(self) -> str:
        """Expression."""
        return ""

    def __getitem__(self, i: int) -> float:
        ...

    def __repr__(self) -> str:
        ...


class VLink:

    """Symbol of links."""

    name: str
    colorSTR: str
    color: Optional[_Color]
    points: Tuple[int, ...]

    def __init__(
        self,
        name: str,
        color_str: str,
        points: tuple,
        color_func: Optional[Callable[..., _Color]] = None
    ):
        ...

    def __contains__(self, point: int) -> bool:
        ...

    def __repr__(self) -> str:
        ...
