# -*- coding: utf-8 -*-

"""Lark parser to parse the expression.

+ PMKS
+ Triangular iteration
"""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import (
    Tuple,
    List,
    Dict,
    Iterator,
    Union,
    Optional,
)
from lark import Lark, Transformer
from lark.lexer import Token
try:
    from .expression import VJoint, VPoint
    from .graph import Graph
except ImportError:
    from expression import VJoint, VPoint
    from graph import Graph

# Color dictionary.
_color_list = {
    'Red': (172, 68, 68),
    'Green': (110, 190, 30),
    'Blue': (68, 120, 172),
    'Cyan': (0, 255, 255),
    'Magenta': (255, 0, 255),
    'Brick-Red': (255, 130, 130),
    'Yellow': (255, 255, 0),
    'Gray': (160, 160, 160),
    'Orange': (225, 165, 0),
    'Pink': (225, 192, 230),
    'Black': (0, 0, 0),
    'White': (255, 255, 255),
    'Dark-Red': (128, 0, 0),
    'Dark-Green': (0, 128, 0),
    'Dark-Blue': (0, 0, 128),
    'Dark-Cyan': (128, 0, 128),
    'Dark-Magenta': (255, 0, 255),
    'Dark-Yellow': (128, 128, 0),
    'Dark-Gray': (128, 128, 128),
    'Dark-Orange': (225, 140, 0),
    'Dark-Pink': (225, 20, 147),
}

color_names = tuple(sorted(_color_list.keys()))


def color_rgb(name: str) -> Tuple[int, int, int]:
    """Get color by name.

    + Invalid color
    + Color key
    + RGB string.
    """
    if not name:
        return 0, 0, 0
    elif name in _color_list:
        return _color_list[name]
    else:
        # Input RGB as a "(255, 255, 255)" string.
        color_text = tuple(int(i) for i in (
            name.replace('(', '')
            .replace(')', '')
            .replace(" ", '')
            .split(',', maxsplit=3)
        ))
        return color_text[:3]


_COLORS = "|".join(f'"{color}"' for color in reversed(color_names))

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
    LETTER: UCASE_LETTER | LCASE_LETTER
    CNAME: ("_" | LETTER) ("_" | LETTER | DIGIT)*

    // White space and new line
    WS: /[ \t\f\r\n]/+
    CR: /\r/
    LF: /\n/
    NEWLINE: (CR? LF)+
    %ignore WS
    %ignore NEWLINE

    // Comment
    COMMENT: "#" /[^\n]/*
    %ignore COMMENT

    // Custom data type
    JOINTTYPE: "RP" | "R" | "P"
    COLOR: """ + _COLORS + """
    type: JOINTTYPE
    name: CNAME
    number: SIGNED_NUMBER
    color_value: INT

    // Main grammar
    joint: "J[" type ["," angle] ["," color] "," point "," link "]"
    link: "L[" name ("," name)* "]"
    point: "P[" number  "," number "]"
    angle: "A[" number "]"
    color: "color[" (("(" color_value ("," color_value) ~ 2 ")") | COLOR) "]"
    mechanism: "M[" [joint ("," joint)* ","?] "]"
""", start='mechanism')


class _PMKSParams(Transformer):

    """Usage: tree = parser.parse(expr)

    args = transformer().transform(tree)
    args: Dict[str, value]
    """

    def type(self, n: List[Token]) -> str:
        return str(n[0])

    name = type

    def color(self, n: List[Token]) -> str:
        return str(n[0]) if len(n) == 1 else str(tuple(n))

    def color_value(self, n: List[Token]) -> int:
        return int(n[0])

    def number(self, n: List[Token]) -> float:
        return float(n[0])

    def point(self, c: List[Token]) -> Tuple[Token]:
        return tuple(c)

    angle = number

    def link(self, a: List[Token]) -> Tuple[Token]:
        return tuple(a)

    def joint(self, args: List[Token]) -> List[Union[str, int, float]]:
        """Sort the argument list.

        [0]: type
        ([1]: angle)
        ([1])[2]: color
        ([2])[3]: point (coordinate)
        ([3])[4]: link
        """
        has_angle = args[0] != 'R'
        x, y = args[-2]
        return [
            ','.join(args[-1]),
            f'{args[0]}:{args[1]}' if has_angle else 'R',
            args[2] if has_angle else args[1],
            x,
            y
        ]

    def mechanism(self, j: List[Token]) -> List[Token]:
        return j


class _PMKSVPoints(_PMKSParams):

    """Using same grammar return as VPoints."""

    def type(self, n: List[Token]) -> VJoint:
        """Return as int type."""
        type_str = str(n[0])
        if type_str == 'R':
            return VJoint.R
        elif type_str == 'P':
            return VJoint.P
        elif type_str == 'RP':
            return VJoint.RP

    def joint(self, args: List[Token]) -> VPoint:
        """Same as parent."""
        has_angle = args[0] != 0
        x, y = args[-2]
        return VPoint(
            ','.join(args[-1]),
            args[0],
            args[1] if has_angle else 0.,
            args[2] if has_angle else args[1],
            x,
            y
        )


def parse_params(expr: str) -> List[List[Union[str, float]]]:
    """Using to parse the expression and return arguments."""
    return _PMKSParams().transform(_GRAMMAR.parse(expr))


def parse_vpoints(expr: str) -> List[VPoint]:
    """Parse as VPoints."""
    return _PMKSVPoints().transform(_GRAMMAR.parse(expr))


def edges_view(graph: Graph) -> Iterator[Tuple[int, Tuple[int, int]]]:
    """This generator can keep the numbering be consistent."""
    yield from enumerate(sorted(tuple(sorted(e)) for e in graph.edges))


def graph2vpoints(
    graph: Graph,
    pos: Dict[int, Tuple[float, float]],
    cus: Optional[Dict[str, int]] = None,
    same: Optional[Dict[int, int]] = None
) -> List[VPoint]:
    """Change NetworkX graph into VPoints.

    cus: custom nodes (not joint)
        {node_name: link_number}
    same: multiple joints
        {n1: n2, n3: n2} => (n1 as n2) and (n3 as n2)
    """
    if cus is None:
        cus: Dict[str, int] = {}
    if same is None:
        same: Dict[int, int] = {}

    same_r = {}
    for k, v in same.items():
        if v in same_r:
            same_r[v].append(k)
        else:
            same_r[v] = [k]
    tmp_list = []
    ev = dict(edges_view(graph))
    for i, e in ev.items():
        if i in same:
            # Do not connect to anyone!
            continue
        e = set(e)
        if i in same_r:
            for j in same_r[i]:
                e.update(set(ev[j]))
        link = ", ".join((str(l) if l else 'ground') for l in e)
        x, y = pos[i]
        tmp_list.append(VPoint.r_joint(link, x, y))
    for name in sorted(cus):
        link = str(cus[name]) if cus[name] else 'ground'
        x, y = pos[int(name.replace('P', ''))]
        tmp_list.append(VPoint.r_joint(link, x, y))
    return tmp_list


try:
    from pygments.lexer import RegexLexer
    from pygments.token import (
        Comment,
        Keyword,
        Name,
        Number,
    )
except ImportError:
    HAS_PYGMENTS = False
else:
    HAS_PYGMENTS = True

    class PMKSLexer(RegexLexer):

        """PMKS highlighter by Pygments."""

        name = 'PMKS'

        tokens = {
            'root': [
                ('#.*$', Comment.Single),
                ('(M)|(J)|(L)|(P)|(A)|(color)', Name.Function),
                ('|'.join(f"({color})" for color in color_names), Name.Variable),
                ('(RP)|(R)|(P)', Keyword.Constant),
                (r'(\d+\.\d*|\d*\.\d+)([eE][+-]?[0-9]+)?j?', Number.Float),
            ]
        }
