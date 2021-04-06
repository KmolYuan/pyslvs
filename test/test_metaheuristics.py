# -*- coding: utf-8 -*-

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from pyslvs.metaheuristics import algorithm, default, AlgorithmType
from pyslvs.metaheuristics.test import TestObj
from . import TestBase


class AlgorithmTest(TestBase):

    def case(self, alg: AlgorithmType):
        """Test with an objective function."""
        settings = {'min_fit': 1e-20, 'report': 10, 'parallel': True}
        obj = TestObj()
        s = default(alg)
        s.update(settings)
        x, fval = algorithm(alg)(obj, s).run()
        self.assertAlmostEqual(0., x[0], 10)
        self.assertAlmostEqual(0., x[1], 10)
        self.assertAlmostEqual(0., fval, 20)

    def test_rga(self):
        self.case(AlgorithmType.RGA)

    def test_de(self):
        self.case(AlgorithmType.DE)

    def test_fa(self):
        self.case(AlgorithmType.FA)

    def test_tlbo(self):
        self.case(AlgorithmType.TLBO)
