# -*- coding: utf-8 -*-

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import TestCase
from timeit import repeat
from pyslvs.metaheuristics import ALGORITHM, PARAMS, AlgorithmType
from pyslvs.metaheuristics.test import TestObj, with_mp, without_mp


class AlgorithmTest(TestCase):

    def test_obj_func(self):
        """Test with an objective function."""
        settings = {'min_fit': 1e-20, 'report': 10, 'parallel': True}
        obj = TestObj()
        for t, setting in PARAMS.items():
            settings.update(setting)
            x, fval = ALGORITHM[t](obj, settings).run()
            self.assertAlmostEqual(0., x[0], 6)
            self.assertAlmostEqual(0., x[1], 6)
            self.assertAlmostEqual(0., fval, 6)


def test_speed(algorithm: AlgorithmType, parallel: bool):
    """Test algorithm performance."""
    settings = {'min_fit': 1e-20, 'report': 10, 'parallel': parallel}
    settings.update(PARAMS[algorithm])
    ALGORITHM[algorithm](TestObj(), settings).run()


if __name__ == '__main__':
    # Time test
    no_mp = min(repeat("f()", number=1, repeat=100, globals={'f': without_mp}))
    mp = min(repeat("f()", number=1, repeat=100, globals={'f': with_mp}))
    print(f"> Radial basis function - is faster? {mp < no_mp}")
    print(f"{no_mp:.020f}s")
    print(f"{mp:.020f}s ... {no_mp / mp * 100 - 100:+.02f}%")
    for alg_type in AlgorithmType:
        no_mp = min(repeat(f"test_speed({alg_type!s}, False)",
                           number=1, repeat=100, globals=globals()))
        mp = min(repeat(f"test_speed({alg_type!s}, True)",
                        number=1, repeat=100, globals=globals()))
        print(f"> {alg_type} - is faster? {mp < no_mp}")
        print(f"{no_mp:.020f}s")
        print(f"{mp:.020f}s ... {no_mp / mp * 100 - 100:+.02f}%")
