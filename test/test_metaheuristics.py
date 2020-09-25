# -*- coding: utf-8 -*-

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import TestCase
from timeit import repeat
from pyslvs.metaheuristics import ALGORITHM, PARAMS, AlgorithmType
from pyslvs.metaheuristics.test import TestObj


class AlgorithmTest(TestCase):

    def test_obj_func(self):
        """Test with an objective function."""
        settings = {'min_fit': 1e-20, 'report': 10}
        obj = TestObj()
        for t, setting in PARAMS.items():
            settings.update(setting)
            x, fval = ALGORITHM[t](obj, settings).run()
            self.assertAlmostEqual(0., x[0], 6)
            self.assertAlmostEqual(0., x[1], 6)
            self.assertAlmostEqual(0., fval, 6)


def test_speed():
    """Test algorithm performance."""
    settings = {'min_fit': 1e-20, 'report': 10}
    obj = TestObj()
    for t, setting in PARAMS.items():
        if t in {AlgorithmType.DE, AlgorithmType.TLBO}:
            continue
        settings.update(setting)
        ALGORITHM[t](obj, settings).run()


if __name__ == '__main__':
    # Time test
    print(min(repeat("test_speed()",
                     "from __main__ import test_speed",
                     number=1, repeat=100)))
