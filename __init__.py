# -*- coding: utf-8 -*-

"""Kernel of Pyslvs.

This kernel can work without GUI.

Modules:
+ Solver:
    + parser
    + tinycadlib
        + Sketch Solve solver
    + triangulation
+ Number synthesis:
    + number
+ Structure Synthesis:
    + atlas
+ Dimensional synthesis:
    + planarlinkage
    + rga
    + firefly
    + de

Dependents:
+ lark-parser
+ pygments (optional: provide highlighting)
"""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2018"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"
__version__ = (18, 12, 0, 'dev')

from .pmks import VPoint, VLink
from .bfgs import vpoint_solving as bfgs_vpoint_solving
from .tinycadlib import (
    Coordinate,
    PLAP,
    PLLP,
    PLPP,
    PXY,
    vpoint_dof,
    expr_parser,
    expr_solving,
    data_collecting,
)
from . import verify
from .planarlinkage import Planar
from .rga import Genetic
from .firefly import Firefly
from .de import Differential
from .number import number_synthesis, contracted_link
from .graph import (
    Graph,
    link_assortments,
    contracted_link_assortments,
)
from .planar_check import is_planar
from .graph_layout import outer_loop_layout
from .atlas import topo
from .triangulation import vpoints_configure
from ._parser import (
    colorNames,
    color_rgb,
    parse_params,
    parse_vpoints,
    HAS_PYGMENTS,
)
from .examples import example_list
if HAS_PYGMENTS:
    from ._parser import PMKSLexer

__all__ = [
    'Genetic',
    'Firefly',
    'Differential',
    'Coordinate',
    'PLAP',
    'PLLP',
    'PLPP',
    'PXY',
    'expr_parser',
    'expr_solving',
    'data_collecting',
    'VPoint',
    'VLink',
    'bfgs_vpoint_solving',
    'Planar',
    'number_synthesis',
    'contracted_link',
    'is_planar',
    'outer_loop_layout',
    'topo',
    'Graph',
    'link_assortments',
    'contracted_link_assortments',
    'vpoints_configure',
    'vpoint_dof',
    'colorNames',
    'color_rgb',
    'parse_params',
    'parse_vpoints',
    'PMKSLexer',
    'example_list',

    # Modules
    'verify',
]
