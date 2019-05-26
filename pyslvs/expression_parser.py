# -*- coding: utf-8 -*-

"""Lark parser to parse the expression."""

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
from lark import Lark, Transformer, LexError
from lark.lexer import Token
try:
    from .expression import (
        get_vlinks,
        VJoint,
        VPoint,
        VLink,
    )
    from .graph import Graph
except ImportError:
    from expression import (
        get_vlinks,
        VJoint,
        VPoint,
        VLink,
    )
    from graph import Graph

# Color dictionary.
_color_list: Dict[str, Tuple[int, int, int]] = {
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
    + RGB string "(r, g, b)".
    """
    if name in _color_list:
        return _color_list[name]
    else:
        try:
            # Input RGB as a "(255, 255, 255)" string.
            color_text: Tuple[int, int, int] = tuple(int(i) for i in (
                name.replace('(', '')
                .replace(')', '')
                .replace(" ", '')
                .split(',', maxsplit=3)
            ))
        except ValueError:
            return 0, 0, 0
        else:
            return color_text


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
    COLOR: """ + "|".join(f'"{color}"' for color in color_names) + r"""
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


class _ParamsTrans(Transformer):

    """Transformer will parse into a list of VPoint data."""

    @staticmethod
    def type(n: List[Token]) -> str:
        return str(n[0])

    name = type

    @staticmethod
    def color(n: List[Token]) -> str:
        return str(n[0]) if len(n) == 1 else str(tuple(n))

    @staticmethod
    def color_value(n: List[Token]) -> int:
        return int(n[0])

    @staticmethod
    def number(n: List[Token]) -> float:
        return float(n[0])

    @staticmethod
    def point(c: List[Token]) -> Tuple[Token]:
        return tuple(c)

    angle = number

    @staticmethod
    def link(a: List[Token]) -> Tuple[Token]:
        return tuple(a)

    @staticmethod
    def joint(args: List[Token]) -> List[Union[str, int, float]]:
        """Sort the argument list.

        [0]: type
        ([1]: angle)
        ([2]: color)
        [-2]: point (coordinate)
        [-1]: link
        """
        type_str = args[0]
        x, y = args[-2]
        links = ','.join(args[-1])
        if type_str == 'R':
            if len(args) == 3:
                return [links, 'R', 'Green', x, y]
            elif len(args) == 4:
                return [links, 'R', args[-3], x, y]
        else:
            type_angle = f'{args[0]}:{args[1]}'
            if len(args) == 4:
                return [links, type_angle, 'Green', x, y]
            elif len(args) == 5:
                return [links, type_angle, args[-3], x, y]

        raise LexError(f"invalid options: {args}")

    @staticmethod
    def mechanism(joints: List[Token]) -> List[Token]:
        return joints


class _PositionTrans(_ParamsTrans):

    """Transformer will parse into a list of position data."""

    @staticmethod
    def joint(args: List[Token]) -> Tuple[float, float]:
        x, y = args[-2]
        return x, y


class _VPointsTrans(_ParamsTrans):

    """Using same grammar return as VPoints."""

    @staticmethod
    def type(n: List[Token]) -> VJoint:
        """Return as int type."""
        type_str = str(n[0])
        if type_str == 'R':
            return VJoint.R
        elif type_str == 'P':
            return VJoint.P
        elif type_str == 'RP':
            return VJoint.RP

    @staticmethod
    def joint(args: List[Token]) -> VPoint:
        """Same as parent."""
        type_int = args[0]
        x, y = args[-2]
        links: Tuple[str, ...] = args[-1]
        if type_int == VJoint.R:
            if len(args) == 3:
                return VPoint.r_joint(links, x, y)
            elif len(args) == 4:
                return VPoint(links, VJoint.R, 0., args[-3], x, y, color_rgb)
        else:
            if len(args) == 4:
                return VPoint.slider_joint(links, type_int, args[1], x, y)
            elif len(args) == 5:
                return VPoint(links, type_int, args[1], args[-3], x, y, color_rgb)

        raise LexError(f"invalid options: {args}")


_params_translator = _ParamsTrans()
_pos_translator = _PositionTrans()
_vpoint_translator = _VPointsTrans()


def parse_params(expr: str) -> List[List[Union[str, float]]]:
    """Using to parse the expression and return arguments."""
    return _params_translator.transform(_GRAMMAR.parse(expr))


def parse_pos(expr: str) -> List[Tuple[float, float]]:
    """Using to parse the expression and return arguments."""
    return _pos_translator.transform(_GRAMMAR.parse(expr))


def parse_vpoints(expr: str) -> List[VPoint]:
    """Parse as VPoints."""
    return _vpoint_translator.transform(_GRAMMAR.parse(expr))


def parse_vlinks(expr: str) -> List[VLink]:
    """Parse as VLinks."""
    return get_vlinks(parse_vpoints(expr))


def edges_view(graph: Graph) -> Iterator[Tuple[int, Tuple[int, int]]]:
    """This generator can keep the numbering be consistent."""
    yield from enumerate(sorted(tuple(sorted(e)) for e in graph.edges))


def graph2vpoints(
    graph: Graph,
    pos: Dict[int, Tuple[float, float]],
    cus: Optional[Dict[int, int]] = None,
    same: Optional[Dict[int, int]] = None,
    grounded: Optional[int] = None
) -> List[VPoint]:
    """Change NetworkX graph into VPoints.

    cus: custom nodes (not joint)
        {node_name: link_number}
    same: multiple joints
        {n1: n2, n3: n2} => (n1 as n2) and (n3 as n2)
    """
    if cus is None:
        cus: Dict[int, int] = {}
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
    for i, edge in ev.items():
        if i in same:
            # Do not connect to anyone!
            continue
        edge = set(edge)
        if i in same_r:
            for j in same_r[i]:
                edge.update(set(ev[j]))
        x, y = pos[i]
        links = [(f"L{link}" if link != grounded else 'ground') for link in edge]
        tmp_list.append(VPoint.r_joint(links, x, y))
    for name in sorted(cus):
        link = f"L{cus[name]}" if cus[name] != grounded else 'ground'
        x, y = pos[name]
        tmp_list.append(VPoint.r_joint((link,), x, y))
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

        """Mechanism expression highlighter by Pygments."""

        name = 'Mechanism Expression'

        tokens = {
            'root': [
                (r'#[^\n]*', Comment.Single),
                (r'#\[[\s\S]*#\][^\n]*', Comment.Multiple),
                (r'M|J|L|P|A|color', Name.Function),
                ('|'.join(f"{color}" for color in reversed(color_names)), Name.Variable),
                (r'RP|R|P', Keyword.Constant),
                (r'(\d+\.\d*|\.\d+)([eE][+-]?\d+)?', Number.Float),
            ]
        }
