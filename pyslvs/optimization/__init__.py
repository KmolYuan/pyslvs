# -*- coding: utf-8 -*-

"""Pyslvs optimization targets."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from .f_planar import (
    FPlanar, norm_path, curvature, derivative, path_signature,
    cross_correlation,
)
from .n_planar import NPlanar, norm_pca
from .utility import FConfig, NConfig
