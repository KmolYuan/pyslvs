# -*- coding: utf-8 -*-

"""Kernel of Pyslvs."""

__version__ = "20.10.0.post0"
__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from .expression import get_vlinks, VJoint, VPoint, VLink, Coord
from .bfgs import SolverSystem
from .topo_config import t_config, EStack
from .tinycadlib import (
    pxy, ppp, plap, pllp, plpp, palp, vpoint_dof, expr_solving,
)
from .planar_linkage import (
    FMatch, norm_path, curvature, derivative, path_signature, cross_correlation,
)
from .expression_parser import (
    color_names, color_rgb, parse_params, parse_pos, parse_vpoints,
    parse_vlinks, edges_view, graph2vpoints, PointArgs, LinkArgs,
)
from .example import example_list, all_examples
from .collection import collection_list, all_collections
from .efd import efd_fitting

# backward alignment
Coordinate = Coord


def get_include() -> str:
    """Get include directory."""
    from os.path import dirname
    return dirname(__file__)
