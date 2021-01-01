# -*- coding: utf-8 -*-

"""Lark parser to parse the expression."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from abc import abstractmethod
from typing import (
    cast, Tuple, List, Dict, Iterator, Optional, Union, TypeVar, Generic,
)
from dataclasses import dataclass
from lark import Lark, Transformer, LexError
from .expression import get_vlinks, VJoint, VPoint, VLink
from .graph import Graph

_T1 = TypeVar('_T1')
_T2 = TypeVar('_T2')
_Coord = Tuple[float, float]
_JointArgs = List[Union[str, VJoint, float, _Coord, Tuple[str, ...]]]

# Color dictionary
_color_list: Dict[str, Tuple[int, int, int]] = {
    'red': (172, 68, 68),
    'green': (110, 190, 30),
    'blue': (68, 120, 172),
    'cyan': (0, 255, 255),
    'magenta': (255, 0, 255),
    'brick-red': (255, 130, 130),
    'yellow': (255, 255, 0),
    'gray': (160, 160, 160),
    'orange': (225, 165, 0),
    'pink': (225, 192, 230),
    'black': (0, 0, 0),
    'white': (255, 255, 255),
    'dark-red': (128, 0, 0),
    'dark-green': (0, 128, 0),
    'dark-blue': (0, 0, 128),
    'dark-cyan': (128, 0, 128),
    'dark-magenta': (255, 0, 255),
    'dark-yellow': (128, 128, 0),
    'dark-gray': (128, 128, 128),
    'dark-orange': (225, 140, 0),
    'dark-pink': (225, 20, 147),
}
color_names = tuple(sorted(_color_list.keys()))


def color_rgb(name: str) -> Tuple[int, int, int]:
    """Get color by name.

    Get RGB color data by name, return `(0, 0, 0)` if it is invalid.
    Also support `"(R, G, B)"` string format.
    """
    name = name.lower()
    if name in _color_list:
        return _color_list[name]
    else:
        try:
            # Input RGB as a "(255, 255, 255)" string
            rgb = (
                name.replace('(', '')
                .replace(')', '')
                .replace(" ", '')
                .split(',', maxsplit=3)
            )
            color_text = (int(rgb[0]), int(rgb[1]), int(rgb[2]))
        except ValueError:
            return 0, 0, 0
        else:
            return color_text


@dataclass(repr=False, eq=False)
class PointArgs:
    """Point table argument."""
    links: str
    type: str
    color: str
    x: float
    y: float


@dataclass(repr=False, eq=False)
class LinkArgs:
    """Link table argument."""
    name: str
    color: str
    points: str


_GRAMMAR = Lark(r"""
// Number
DIGIT: "0".."9"
INT: DIGIT+
SIGNED_INT: ["+" | "-"] INT
DECIMAL: INT "." INT? | "." INT
_EXP: ("e" | "E") SIGNED_INT
FLOAT: INT _EXP | DECIMAL _EXP?
NUMBER: FLOAT | INT
SIGNED_NUMBER: ["+" | "-"] NUMBER

// Letters
LCASE_LETTER: "a".."z"
UCASE_LETTER: "A".."Z"
LETTER: UCASE_LETTER | LCASE_LETTER | "_"
CNAME: LETTER (LETTER | DIGIT)*

// White space and new line
WS: /\s+/
CR: /\r/
LF: /\n/
NEWLINE: (CR? LF)+
%ignore WS
%ignore NEWLINE

// Comment
LINE_COMMENT: /#[^\n]*/
MULTILINE_COMMENT: /#\[[\s\S]*#\][^\n]*/
%ignore LINE_COMMENT
%ignore MULTILINE_COMMENT

// Custom data type
JOINT_TYPE: "RP" | "R" | "P"
COLOR: """ + "|".join(f'"{color}"i' for color in color_names) + r"""
type: JOINT_TYPE
name: CNAME
number: SIGNED_NUMBER
color_value: INT

// Main grammar
joint: "J[" type ["," angle] ["," color] "," point "," link "]"
link: "L[" [name ("," name)* ","?] "]"
point: "P[" number "," number "]"
angle: "A[" number "]"
color: "color[" (("(" color_value ("," color_value) ~ 2 ")") | COLOR) "]"
mechanism: "M[" [joint ("," joint)* ","?] "]"
?start: mechanism
""", parser='lalr')


class _Transformer(Transformer, Generic[_T1, _T2]):
    """Base transformer implementation."""

    @staticmethod
    @abstractmethod
    def type(n: List[str]) -> _T1:
        raise NotImplementedError

    @staticmethod
    def name(n: List[str]) -> str:
        return str(n[0])

    @staticmethod
    def color(n: List[str]) -> str:
        return str(n[0]) if len(n) == 1 else str(tuple(n))

    @staticmethod
    def color_value(n: List[str]) -> int:
        return int(n[0])

    @staticmethod
    def number(n: List[str]) -> float:
        return float(n[0])

    angle = number

    @staticmethod
    def point(c: List[float]) -> _Coord:
        return c[0], c[1]

    @staticmethod
    def link(a: List[str]) -> Tuple[str, ...]:
        return tuple(a)

    @staticmethod
    @abstractmethod
    def joint(args: _JointArgs) -> _T2:
        raise NotImplementedError

    @staticmethod
    def mechanism(joints: List[_T2]) -> List[_T2]:
        return joints


class _ParamsTrans(_Transformer[str, PointArgs]):
    """Transformer will parse into a list of VPoint data."""

    @staticmethod
    def type(n: List[str]) -> str:
        return str(n[0])

    @staticmethod
    def joint(args: _JointArgs) -> PointArgs:
        """Sort the argument list.

        [0]: type
        ([1]: angle)
        ([2]: color)
        [-2]: point (coordinate)
        [-1]: link
        """
        type_str = cast(str, args[0])
        x, y = cast(_Coord, args[-2])
        links = ','.join(cast(Tuple[str, ...], args[-1]))
        if type_str == 'R':
            if len(args) == 3:
                return PointArgs(links, 'R', 'Green', x, y)
            elif len(args) == 4:
                color = cast(str, args[-3])
                return PointArgs(links, 'R', color, x, y)
        else:
            angle = cast(float, args[1])
            type_angle = f'{type_str}:{angle}'
            if len(args) == 4:
                return PointArgs(links, type_angle, 'Green', x, y)
            elif len(args) == 5:
                color = cast(str, args[-3])
                return PointArgs(links, type_angle, color, x, y)

        raise LexError(f"invalid options: {args}")


class _PositionTrans(_Transformer[str, _Coord]):
    """Transformer will parse into a list of position data."""

    @staticmethod
    def type(n: List[str]) -> str:
        return str(n[0])

    @staticmethod
    def joint(args: _JointArgs) -> _Coord:
        x, y = cast(_Coord, args[-2])
        return x, y


class _VPointsTrans(_Transformer[VJoint, VPoint]):
    """Using same grammar return as VPoints."""

    @staticmethod
    def type(n: List[str]) -> VJoint:
        """Return as int type."""
        type_str = str(n[0])
        if type_str == 'R':
            return VJoint.R
        elif type_str == 'P':
            return VJoint.P
        elif type_str == 'RP':
            return VJoint.RP
        else:
            raise ValueError(f"invalid joint type: {type_str}")

    @staticmethod
    def joint(args: _JointArgs) -> VPoint:
        """Sort the argument list.

        [0]: type
        ([1]: angle)
        ([2]: color)
        [-2]: point (coordinate)
        [-1]: link
        """
        type_int = cast(VJoint, args[0])
        x, y = cast(_Coord, args[-2])
        links = cast(Tuple[str, ...], args[-1])
        if type_int == VJoint.R:
            if len(args) == 3:
                return VPoint.r_joint(links, x, y)
            elif len(args) == 4:
                color = cast(str, args[-3])
                return VPoint(links, VJoint.R, 0., color, x, y, color_rgb)
        else:
            angle = cast(float, args[1])
            if len(args) == 4:
                return VPoint.slider_joint(links, type_int, angle, x, y)
            elif len(args) == 5:
                color = cast(str, args[-3])
                return VPoint(links, type_int, angle, color, x, y, color_rgb)

        raise LexError(f"invalid options: {args}")


_params_translator = _ParamsTrans()
_pos_translator = _PositionTrans()
_vpoint_translator = _VPointsTrans()


def parse_params(expr: str) -> List[PointArgs]:
    """Parse mechanism expression into VPoint constructor arguments."""
    return _params_translator.transform(_GRAMMAR.parse(expr))


def parse_pos(expr: str) -> List[_Coord]:
    """Parse mechanism expression into coordinates."""
    return _pos_translator.transform(_GRAMMAR.parse(expr))


def parse_vpoints(expr: str) -> List[VPoint]:
    """Parse mechanism expression into VPoint objects."""
    return _vpoint_translator.transform(_GRAMMAR.parse(expr))


def parse_vlinks(expr: str) -> List[VLink]:
    """Parse mechanism expression into VLink objects."""
    return get_vlinks(parse_vpoints(expr))


def _sorted_pair(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def edges_view(graph: Graph) -> Iterator[Tuple[int, Tuple[int, int]]]:
    """The iterator will yield the sorted edges from `graph`."""
    yield from enumerate(sorted(_sorted_pair(n1, n2) for n1, n2 in graph.edges))


def graph2vpoints(
    graph: Graph,
    pos: Dict[int, _Coord],
    cus: Optional[Dict[int, int]] = None,
    same: Optional[Dict[int, int]] = None,
    grounded: Optional[int] = None
) -> List[VPoint]:
    """Transform `graph` into [VPoint] objects. The vertices are mapped to links.

    + `pos`: Position for each vertices.
    + `cus`: Extra points on the specific links.
    + `same`: Multiple joint setting. The joints are according to [`edges_view`](#edges_view).
    + `grounded`: The ground link of vertices.
    """
    if cus is None:
        cus = {}
    if same is None:
        same = {}
    same_r: Dict[int, List[int]] = {}
    for k, v in same.items():
        if v in same_r:
            same_r[v].append(k)
        else:
            same_r[v] = [k]
    tmp_list = []
    ev = dict(edges_view(graph))
    for i, edge in ev.items():
        if i in same:
            # Do not connect to anyone!
            continue
        edges = set(edge)
        if i in same_r:
            for j in same_r[i]:
                edges.update(set(ev[j]))
        x, y = pos[i]
        links = [
            f"L{link}" if link != grounded else VLink.FRAME for link in edges
        ]
        tmp_list.append(VPoint.r_joint(links, x, y))
    for name in sorted(cus):
        link = f"L{cus[name]}" if cus[name] != grounded else VLink.FRAME
        x, y = pos[name]
        tmp_list.append(VPoint.r_joint((link,), x, y))
    return tmp_list
