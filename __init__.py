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
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"
__version__ = (19, 5, 0, 'dev')
__version_str__ = f"{__version__[0]}.{__version__[1]:02d}.{__version__[2]} ({__version__[3]})"

from .expression import (
    get_vlinks,
    VJoint,
    VPoint,
    VLink,
    Coordinate,
)
from .bfgs import vpoint_solving
from .triangulation import vpoints_configure, ExpressionStack
from .tinycadlib import (
    plap,
    pllp,
    plpp,
    pxy,
    vpoint_dof,
    expr_parser,
    expr_solving,
    data_collecting,
)
from .verify import Verification
from .planar_linkage import Planar
from .rga import Genetic
from .firefly import Firefly
from .de import Differential
from .number import link_synthesis, contracted_link_synthesis
from .graph import (
    Graph,
    link_assortment,
    contracted_link_assortment,
    labeled_enumerate,
)
from .planar_check import is_planar
from .graph_layout import external_loop_layout
from .atlas import conventional_graph, contracted_graph
from .expression_parser import (
    color_names,
    color_rgb,
    parse_params,
    parse_pos,
    parse_vpoints,
    parse_vlinks,
    edges_view,
    graph2vpoints,
    HAS_PYGMENTS,
)
from .example import example_list
from .collection import collection_list
if HAS_PYGMENTS:
    from .expression_parser import PMKSLexer

__all__ = [
    '__version__',
    '__version_str__',
    'Genetic',
    'Firefly',
    'Differential',
    'Coordinate',
    'plap',
    'pllp',
    'plpp',
    'pxy',
    'expr_parser',
    'expr_solving',
    'data_collecting',
    'get_vlinks',
    'VJoint',
    'VPoint',
    'VLink',
    'vpoint_solving',
    'Planar',
    'link_synthesis',
    'contracted_link_synthesis',
    'is_planar',
    'external_loop_layout',
    'conventional_graph',
    'contracted_graph',
    'Graph',
    'link_assortment',
    'contracted_link_assortment',
    'labeled_enumerate',
    'vpoints_configure',
    'ExpressionStack',
    'vpoint_dof',
    'color_names',
    'color_rgb',
    'parse_params',
    'parse_pos',
    'parse_vpoints',
    'parse_vlinks',
    'edges_view',
    'graph2vpoints',
    'PMKSLexer',
    'example_list',
    'collection_list',
    'Verification',
]
