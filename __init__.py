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
__version__ = (19, 1, 0, 'dev')

import numpy.core._dtype_ctypes
from .expression import VJoint, VPoint, VLink
from .bfgs import vpoint_solving
from .tinycadlib import (
    Coordinate,
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
from .number import number_synthesis, contracted_link
from .graph import (
    Graph,
    link_assortments,
    contracted_link_assortments,
)
from .planar_check import is_planar
from .graph_layout import external_loop_layout
from .atlas import topo
from .triangulation import vpoints_configure
from ._parser import (
    color_names,
    color_rgb,
    parse_params,
    parse_vpoints,
    edges_view,
    graph2vpoints,
    HAS_PYGMENTS,
)
from .examples import example_list
if HAS_PYGMENTS:
    from ._parser import PMKSLexer

__all__ = [
    '__version__',
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
    'VJoint',
    'VPoint',
    'VLink',
    'vpoint_solving',
    'Planar',
    'number_synthesis',
    'contracted_link',
    'is_planar',
    'external_loop_layout',
    'topo',
    'Graph',
    'link_assortments',
    'contracted_link_assortments',
    'vpoints_configure',
    'vpoint_dof',
    'color_names',
    'color_rgb',
    'parse_params',
    'parse_vpoints',
    'edges_view',
    'graph2vpoints',
    'PMKSLexer',
    'example_list',
    'Verification',
]
