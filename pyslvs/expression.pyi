# -*- coding: utf-8 -*-

from typing import Tuple, List, Iterable, Sequence, Callable, Optional, ClassVar
from enum import IntEnum, auto

# Color type must be a RGB data
_Color = Tuple[int, int, int]
_Coord = Tuple[float, float]

def get_vlinks(vpoints: Iterable[VPoint]) -> List[VLink]:
    """Get VLinks of a list of VPoint."""
    ...

class Coordinate:

    """A class to store the coordinate."""

    x: float
    y: float

    def __init__(self, x: float, y: float) -> None:
        ...

    def distance(self, p: Coordinate) -> float:
        """Distance."""
        ...

    def is_nan(self) -> bool:
        """Test this coordinate is a error-occurred answer."""
        ...

    def __repr__(self) -> str:
        ...

class VJoint(IntEnum):
    R = auto()
    P = auto()
    RP = auto()

class VPoint:

    """Symbol of joints."""

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
        ...

    @staticmethod
    def r_joint(links: Iterable[str], x: float, y: float) -> VPoint:
        """Create by coordinate."""
        ...

    @staticmethod
    def slider_joint(links: Iterable[str], type_int: VJoint, angle: float, x: float, y: float) -> VPoint:
        """Create by coordinate."""
        ...

    def copy(self) -> VPoint:
        """Copy method of Python."""
        ...

    @property
    def cx(self) -> float:
        """X value of first current coordinate."""
        return 0.

    @property
    def cy(self) -> float:
        """Y value of first current coordinate."""
        return 0.

    def set_links(self, links: Iterable[str]) -> None:
        """Set links."""
        ...

    def replace_link(self, link1: str, link2: str) -> None:
        """Replace link1 as link2."""
        ...

    def move(self, c1: Tuple[float, float], c2: Optional[Tuple[float, float]] = None) -> None:
        """Change coordinates of this point."""
        ...

    def locate(self, x: float, y: float) -> None:
        """Change the origin coordinate of this point directly."""
        ...

    def rotate(self, angle: float) -> None:
        """Change the angle of slider slot by degrees."""
        ...

    def set_offset(self, offset: float) -> None:
        """Set slider offset."""
        ...

    def disable_offset(self) -> None:
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

    def expr(self) -> str:
        """Expression."""
        ...

    def __getitem__(self, i: int) -> float:
        ...

    def __repr__(self) -> str:
        ...

class VLink:

    """Symbol of links."""

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
        ...

    def set_points(self, points: Iterable[int]) -> None:
        """Set points."""
        ...

    def __contains__(self, point: int) -> bool:
        ...

    def __repr__(self) -> str:
        ...
