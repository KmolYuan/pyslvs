# -*- coding: utf-8 -*-

"""Pyslvs planar linkage module test."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import TestCase
from math import radians, hypot, sin, cos
from numpy import array
from pyslvs import parse_vpoints, norm_path, collection_list, FMatch
from pyslvs.metaheuristics import ALGORITHM, AlgorithmType, PARAMS

_FOUR_BAR = collection_list("Four bar linkage mechanism")
_FOUR_BAR.update({
    'expression': parse_vpoints(_FOUR_BAR['expression']),
    'placement': {0: (-70, -70, 10), 1: (70, -70, 10)},
    'target': {
        4: [
            (60.3, 118.12),
            (31.02, 115.62),
            (3.52, 110.62),
            (-25.77, 104.91),
            (-81.49, 69.19),
            (-96.47, 54.906),
            (-109.34, 35.98),
            (-121.84, 13.83),
            (-127.56, -20.09),
            (-128.63, -49.74),
            (-117.56, -65.45),
        ]
    },
    'upper': 100.,
    'lower': 0.,
})
PLANAR_OBJECT = FMatch(_FOUR_BAR)
PATH = [
    (6.7700907146387586, 24.644877369732),
    (3.9327689792658944, 26.12795413801081),
    (-0.907462683602656, 26.05153570209126),
    (-6.24548864186997, 24.8225463940994),
    (-11.133346308106713, 23.18479875363321),
    (-15.520041274243832, 21.57064857277561),
    (-19.801999741184225, 19.85946369468313),
    (-24.175969350603772, 17.70314933508058),
    (-28.438413573857005, 14.9864798961492),
    (-32.292254757096046, 11.90810680124495),
    (-35.67126221930281, 8.6512747009891),
    (-38.64892957450998, 5.10259170466386),
    (-41.047776013814655, 1.021821684545207),
    (-42.2630865594439, -3.473263961136975),
    (-41.608226040938476, -7.65939066663121),
    (-38.904742969958264, -10.5580841883062),
    (-34.73026075172416, -11.62017315187812),
    (-30.055836486064223, -11.09019975630983),
    (-25.618933599062636, -9.68103804133981),
    (-21.627062086461947, -7.87802719739308),
    (-18.00286942796385, -5.56985492705963),
    (-14.787669532074853, -2.3636449904643),
    (-12.175023083437713, 1.772565914254553),
    (-10.11835329911016, 6.2409254685245),
    (-8.025279517469247, 10.22777039213621),
    (-5.057922208547721, 13.38647752586898),
    (-0.9343218889706311, 16.0769324633211),
    (3.4945236042052983, 18.8944943173922),
    (6.54573633651432, 21.9527700520141),
    (6.770090714638762, 24.64487736973219),
]


class PlanarTest(TestCase):

    def test_path_normalization(self):
        """Test path normalization function."""
        path1 = norm_path(PATH)
        alpha = radians(50)
        c = cos(alpha)
        s = sin(alpha)
        path2 = norm_path(array(path1) @ array([[c, -s], [s, c]]))
        for i in range(len(PATH)):
            h = hypot(path1[i][0] - path2[i][0], path1[i][1] - path2[i][1])
            self.assertAlmostEqual(0, h, 6)

    def algorithm_generic(self, t: AlgorithmType):
        """Generic algorithm setup."""
        settings = {'max_gen': 10, 'report': 10}
        settings.update(PARAMS[t])
        algorithm = ALGORITHM[t](PLANAR_OBJECT, settings)
        algorithm.run()
        t_f = algorithm.history()
        self.assertEqual(10, t_f[1][0] - t_f[0][0])

    def test_algorithms(self):
        """Test algorithms."""
        self.assertFalse(PLANAR_OBJECT.is_two_kernel())
        # Real-coded genetic algorithm
        self.algorithm_generic(AlgorithmType.RGA)
        # Firefly algorithm
        self.algorithm_generic(AlgorithmType.Firefly)
        # Differential evolution
        self.algorithm_generic(AlgorithmType.DE)
        # Teaching learning based optimization
        self.algorithm_generic(AlgorithmType.TLBO)
