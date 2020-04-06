# -*- coding: utf-8 -*-

"""Kernel of Pyslvs."""

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
    'norm_path',
    'Planar',
    't_config',
    'EStack',
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
    'all_examples',
    'collection_list',
    'all_collections',
    'efd_fitting',
    'get_include',
]
__version__ = "20.03.0.post0"
__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

import pywt as _
from .expression import get_vlinks, VJoint, VPoint, VLink, Coordinate
from .bfgs import SolverSystem
from .triangulation import t_config, EStack
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
from .planar_linkage import Planar, norm_path
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
from .example import example_list, all_examples
from .collection import collection_list, all_collections
from .efd import efd_fitting


def get_include() -> str:
    """Get include directory."""
    from os.path import dirname
    return dirname(__file__)
