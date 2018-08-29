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
+ Topologic synthesis:
    + topologic
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
__version__ = (18, 9, 0, 'dev')

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
from .verify import Verification
from .planarlinkage import Planar
from .rga import Genetic
from .firefly import Firefly
from .de import DiffertialEvolution
from .number import number_synthesis
from .topologic import topo, Graph
from .triangulation import vpoints_configure
from ._parser import (
    colorNames,
    colorRGB,
    parse_params,
    parse_vpoints,
    HAS_PYGMENTS,
)
if HAS_PYGMENTS:
    from ._parser import PMKSLexer
from .examples import example_list

__all__ = [
    'Genetic',
    'Firefly',
    'DiffertialEvolution',
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
    'Verification',
    'Planar',
    'number_synthesis',
    'topo',
    'Graph',
    'vpoints_configure',
    'vpoint_dof',
    'colorNames',
    'colorRGB',
    'parse_params',
    'parse_vpoints',
    'PMKSLexer',
    'example_list',
]
