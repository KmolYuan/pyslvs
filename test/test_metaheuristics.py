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
        settings = {'min_fit': 1e-20, 'report': 10}
        obj = TestObj()
        s = default(alg)
        s.update(settings)
        a = algorithm(alg)(obj, s)
        ans = a.run()
        x, y = a.result()
        self.assertTrue(ans < 1e-20)
        for i in range(4):
            self.assertTrue(x[i] < 1e-10)
        self.assertEqual(y, ans)

    def test_rga(self):
        self.case(AlgorithmType.RGA)

    def test_de(self):
        self.case(AlgorithmType.DE)

    def test_pso(self):
        self.case(AlgorithmType.PSO)

    def test_fa(self):
        self.case(AlgorithmType.FA)

    def test_tlbo(self):
        self.case(AlgorithmType.TLBO)
