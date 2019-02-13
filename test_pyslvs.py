# -*- coding: utf-8 -*-

"""This module will test the functions of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

import unittest
from unittest import TestCase
from math import sqrt, radians
from bfgs import vpoint_solving
from tinycadlib import (
    Coordinate,
    plap,
    pllp,
    plpp,
    pxy,
    expr_solving,
    data_collecting,
)
from planar_linkage import Planar
from rga import Genetic
from firefly import Firefly
from de import Differential
from number import number_synthesis, contracted_link
from atlas import topo
from graph import (
    Graph,
    link_assortments,
    contracted_link_assortments,
)
from planar_check import is_planar
from graph_layout import external_loop_layout
from triangulation import vpoints_configure
from _parser import parse_vpoints
from example import example_list


_planar_object = Planar({
    'Driver': {'P0': (-70, -70, 50)},
    'Follower': {'P1': (70, -70, 50)},
    'Target': {
        'P4': [
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
    'Expression':
        "plap[P0,L0,a0](P2);"
        "pllp[P2,L1,L2,P1](P3);"
        "pllp[P2,L3,L4,P3](P4)",
    'upper': [100., 100., 100., 100., 100., 360.],
    'lower': [5., 5., 5., 5., 5., 0.],
})


class CoreTest(TestCase):
    """Testing Cython libs."""

    def test_plap(self):
        """Test for plap function."""
        x, y = plap(Coordinate(0, 0), 50 * sqrt(2), radians(45), Coordinate(50, 0))
        self.assertAlmostEqual(50, x)
        self.assertAlmostEqual(50, y)

    def test_pllp(self):
        """Test for pllp function."""
        c1 = Coordinate(-30, 0)
        c2 = Coordinate(30, 0)
        x, y = pllp(c1, 50, 50, c2)
        self.assertAlmostEqual(0, x)
        self.assertAlmostEqual(40, y)
        x, y = pllp(c1, 30, 30, c2)
        self.assertAlmostEqual(x, 0)
        self.assertAlmostEqual(y, 0)
        x, y = pllp(c1, 90, 30, c2)
        self.assertAlmostEqual(60, x)
        self.assertAlmostEqual(0, y)

    def test_plpp(self):
        """Test for plpp function."""
        x, y = plpp(Coordinate(0, 0), sqrt(5), Coordinate(0, -3), Coordinate(3 / 2, 0))
        self.assertAlmostEqual(2, x)
        self.assertAlmostEqual(1, y)

    def test_pxy(self):
        x, y = pxy(Coordinate(80, 90), 40, -20)
        self.assertAlmostEqual(120, x)
        self.assertAlmostEqual(70, y)

    def test_graph_function(self):
        """Test 'graph' libraries."""
        g1 = Graph([(0, 1), (0, 4), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)])
        self.assertFalse(g1.is_degenerate())
        self.assertTrue(g1.is_connected())
        self.assertEqual(1, g1.dof())
        g2 = Graph([(0, 2), (0, 4), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)])
        g3 = Graph([(0, 1), (0, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5)])
        self.assertTrue(g1.is_isomorphic(g2))
        self.assertFalse(g1.is_isomorphic(g3))

        g1 = Graph([(0, 1), (2, 7), (1, 5), (1, 6), (3, 6), (0, 4), (3, 7), (2, 5), (3, 4), (0, 2)])
        self.assertTrue(is_planar(g1))
        self.assertEqual([4, 4], link_assortments(g1))
        self.assertEqual([4, 0, 0, 0], contracted_link_assortments(g1))

        g1 = Graph([(0, 1), (4, 10), (1, 3), (2, 9), (5, 6), (4, 5), (5, 7), (8, 10),
                    (1, 8), (9, 11), (3, 6), (0, 4), (3, 7), (2, 5), (0, 2), (4, 11)])
        pos = external_loop_layout(g1, True)
        self.assertEqual(set(g1.nodes), set(pos))

    def test_atlas(self):
        """Test 'atlas' libraries."""
        answers = []

        type_1 = [4, 2]
        for c_j in contracted_link(type_1):
            answer, _ = topo(type_1, c_j)
            answers.extend(answer)
        self.assertEqual(2, len(answers))
        answers.clear()

        answers_degenerated = []
        for type_2 in ([4, 4, 0], [5, 2, 1], [6, 0, 2]):
            for c_j in contracted_link(type_2):
                answer, _ = topo(type_2, c_j)
                answers.extend(answer)
                answer, _ = topo(type_2, c_j, 2)
                answers_degenerated.extend(answer)
        self.assertEqual(16, len(answers))
        self.assertEqual(40, len(answers_degenerated))

    def test_solving(self):
        """Test triangular formula solving.

        + Test for PMKS parser.
        + Test data collecting function.
        + Test expression solving function.
        """
        expr, inputs = example_list["Jansen's linkage (Single)"]
        vpoints = parse_vpoints(expr)
        self.assertTrue(len(vpoints) == 8)
        exprs = vpoints_configure(vpoints, inputs)
        mapping = {n: f'P{n}' for n in range(len(vpoints))}
        data_dict, dof = data_collecting(exprs, mapping, vpoints)
        for link, link_length in (
            ('L0', 15.002083),
            ('L1', 41.501876),
            ('L2', 49.994906),
            ('L3', 40.09651),
            ('L4', 55.802532),
            ('L5', 61.905252),
            ('L6', 39.3028),
            ('L7', 36.697676),
            ('L8', 39.395233),
            ('L9', 48.995887),
            ('L10', 65.699401),
        ):
            self.assertAlmostEqual(data_dict[link], link_length, 6)
        self.assertEqual(1, dof)
        result = expr_solving(exprs, mapping, vpoints, [0.])
        x, y = result[-1]
        self.assertAlmostEqual(-43.170055, x, 6)
        self.assertAlmostEqual(-91.753226, y, 6)

    def test_solving_bfgs(self):
        """Test Sketch Solve kernel."""
        expr, _ = example_list["Jansen's linkage (Single)"]
        vpoints = parse_vpoints(expr)
        result = vpoint_solving(vpoints, {(0, 1): 0.})
        x, y = result[-1]
        self.assertAlmostEqual(-43.170055, x, 6)
        self.assertAlmostEqual(-91.753226, y, 6)

    def test_number_synthesis(self):
        """Test Number Synthesis function."""
        for nl, nj in [(4, 4), (6, 7), (8, 9), (10, 12)]:
            for factors in number_synthesis(nl, nj):
                count = 0
                for i, factor in enumerate(factors):
                    count += factor * (i + 2)
                self.assertEqual(nj, int(count / 2))

    def test_algorithm_rga(self):
        """Real-coded genetic algorithm."""
        _, t_f = Genetic(_planar_object, {
            'max_time': 1,
            'report': 10,
            # Genetic
            'nPop': 500,
            'pCross': 0.95,
            'pMute': 0.05,
            'pWin': 0.95,
            'bDelta': 5.,
        }).run()
        self.assertTrue(10 >= t_f[1][0] - t_f[0][0])

    def test_algorithm_firefly(self):
        """Firefly algorithm."""
        _, t_f = Firefly(_planar_object, {
            'max_time': 1,
            'report': 10,
            # Firefly
            'n': 80,
            'alpha': 0.01,
            'betaMin': 0.2,
            'gamma': 1.,
            'beta0': 1.,
        }).run()
        self.assertTrue(10 >= t_f[1][0] - t_f[0][0])

    def test_algorithm_de(self):
        """Differential evolution."""
        _, t_f = Differential(_planar_object, {
            'max_time': 1,
            'report': 10,
            # DE
            'strategy': 1,
            'NP': 400,
            'F': 0.6,
            'CR': 0.9,
        }).run()
        self.assertTrue(10 >= t_f[1][0] - t_f[0][0])


if __name__ == '__main__':
    unittest.main()
