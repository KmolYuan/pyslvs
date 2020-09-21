# -*- coding: utf-8 -*-

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import TestCase
from pyslvs.metaheuristics import ALGORITHM, PARAMS, AlgorithmType
from pyslvs.metaheuristics.test import TestObj


class AlgorithmTest(TestCase):

    def test_obj_func(self):
        """Test with an objective function."""
        settings = {'min_fit': 1e-20, 'report': 10}
        obj = TestObj()
        for t, setting in PARAMS.items():
            if t == AlgorithmType.RGA:
                continue
            settings.update(setting)
            x, fval = ALGORITHM[t](obj, settings).run()
            self.assertAlmostEqual(0., x[0], 6)
            self.assertAlmostEqual(0., x[1], 6)
            self.assertAlmostEqual(0., fval, 6)
