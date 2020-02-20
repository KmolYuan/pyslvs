# -*- coding: utf-8 -*-

"""This module contains test objects."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from pyslvs import collection_list, parse_vpoints, Planar

_FOUR_BAR = collection_list("Four bar linkage mechanism")
_FOUR_BAR.update({
    'expression': parse_vpoints(_FOUR_BAR['expression']),
    'placement': {0: (-70, -70, 50), 1: (70, -70, 50)},
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
    'upper': [100., 100., 100., 100., 100.],
    'lower': [0., 0., 0., 0., 0.],
})
PLANAR_OBJECT = Planar(_FOUR_BAR)
DEGREE_CODE_TABLE = [
    (7, [(0, 1), (1, 2), (0, 2)]),
    (51, [(0, 1), (1, 2), (2, 3), (0, 3)]),
    (62, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2)]),
    (63, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 2), (1, 3)]),
    (787, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]),
    (937, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 2)]),
    (504, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (4, 2)]),
    (947, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (4, 2), (1, 4)]),
    (1010, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 2), (0, 3)]),
    (1016, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (4, 2), (0, 2)]),
    (1011, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 2), (0, 3), (1, 4)]),
    (1020, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (4, 2), (0, 2), (1, 4)]),
    (1022,
     [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (4, 2), (0, 2), (1, 4), (1, 3)]),
    (24851, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]),
    (15169, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 5), (3, 5)]),
    (27050, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 2), (3, 5)]),
    (29459, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 2), (1, 3)]),
    (29326, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 3), (1, 5)]),
    (31497, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 3), (0, 4)]),
    (31064, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 2), (0, 4)]),
    (24344, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 3), (0, 5), (2, 5)]),
    (31553, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (0, 2), (0, 5), (2, 5)]),
    (16320, [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (2, 4), (0, 5), (2, 5)]),
    (29327, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 3), (1, 5),
             (2, 4)]),
    (30358, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 4), (1, 5),
             (2, 4)]),
    (31507, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 4), (1, 5),
             (3, 5)]),
    (30485, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (0, 3), (2, 5),
             (2, 4)]),
]
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
