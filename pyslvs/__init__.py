# -*- coding: utf-8 -*-

"""Kernel of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"
__version__ = "20.02.0.dev0"

from .expression import get_vlinks, VJoint, VPoint, VLink, Coordinate
from .bfgs import SolverSystem
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
from . import metaheuristics as _  # Preload module
from .planar_linkage import Planar
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
    PointArgs,
    LinkArgs,
)
from .example import example_list
from .collection import collection_list
from .efd import efd_fitting

__all__ = [
    '__version__',
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
    'SolverSystem',
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
    'PointArgs',
    'LinkArgs',
    'example_list',
    'collection_list',
    'efd_fitting',
]
