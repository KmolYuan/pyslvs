# -*- coding: utf-8 -*-

"""This module will test the functions of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2019"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

import unittest
from unittest import TestCase
from math import sqrt, radians
from copy import deepcopy
from pyslvs import (
    Coordinate,
    vpoint_solving,
    plap,
    pllp,
    plpp,
    pxy,
    expr_solving,
    data_collecting,
    Planar,
    Genetic,
    Firefly,
    Differential,
    link_synthesis,
    contracted_link_synthesis,
    conventional_graph,
    contracted_graph,
    Graph,
    link_assortment,
    contracted_link_assortment,
    is_planar,
    external_loop_layout,
    vpoints_configure,
    parse_vpoints,
    example_list,
    collection_list,
)


_four_bar = deepcopy(collection_list["Four bar linkage mechanism"])
_four_bar['Expression'] = parse_vpoints(_four_bar['Expression'])
_four_bar.update({
    'Placement': {
        0: (-70, -70, 50),
        1: (70, -70, 50),
    },
    'Target': {
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
    'upper': [100., 100., 100., 100., 100., 360.],
    'lower': [0., 0., 0., 0., 0., 0.],
})
_planar_object = Planar(_four_bar)


class CoreTest(TestCase):

    def test_plap(self):
        """Test for plap function."""
        coord = plap(Coordinate(0, 0), 50 * sqrt(2), radians(45), Coordinate(50, 0))
        self.assertAlmostEqual(50, coord.x)
        self.assertAlmostEqual(50, coord.y)

    def test_pllp(self):
        """Test for pllp function."""
        c1 = Coordinate(-30, 0)
        c2 = Coordinate(30, 0)
        coord = pllp(c1, 50, 50, c2)
        self.assertAlmostEqual(0, coord.x)
        self.assertAlmostEqual(40, coord.y)
        coord = pllp(c1, 30, 30, c2)
        self.assertAlmostEqual(coord.x, 0)
        self.assertAlmostEqual(coord.y, 0)
        coord = pllp(c1, 90, 30, c2)
        self.assertAlmostEqual(60, coord.x)
        self.assertAlmostEqual(0, coord.y)

    def test_plpp(self):
        """Test for plpp function."""
        coord = plpp(Coordinate(0, 0), sqrt(5), Coordinate(0, -3), Coordinate(3 / 2, 0))
        self.assertAlmostEqual(2, coord.x)
        self.assertAlmostEqual(1, coord.y)

    def test_pxy(self):
        coord = pxy(Coordinate(80, 90), 40, -20)
        self.assertAlmostEqual(120, coord.x)
        self.assertAlmostEqual(70, coord.y)

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
        self.assertEqual([4, 4], link_assortment(g1))
        self.assertEqual([4, 0, 0, 0], contracted_link_assortment(g1))

        g1 = Graph([(0, 1), (4, 10), (1, 3), (2, 9), (5, 6), (4, 5), (5, 7), (8, 10),
                    (1, 8), (9, 11), (3, 6), (0, 4), (3, 7), (2, 5), (0, 2), (4, 11)])
        pos = external_loop_layout(g1, True)
        self.assertEqual(set(g1.nodes), set(pos))

    def test_atlas(self):
        """Test 'atlas' libraries."""
        answers = []

        type_0 = [4]
        cg_list = contracted_graph(type_0)
        for c_j in contracted_link_synthesis(type_0):
            answer = conventional_graph(cg_list, c_j)
            answers.extend(answer)
        self.assertEqual(1, len(answers))
        answers.clear()

        type_1 = [4, 2]
        cg_list = contracted_graph(type_1)
        for c_j in contracted_link_synthesis(type_1):
            answer = conventional_graph(cg_list, c_j)
            answers.extend(answer)
        self.assertEqual(2, len(answers))
        answers.clear()

        answers_degenerated = []
        for type_2 in ([4, 4, 0], [5, 2, 1], [6, 0, 2]):
            cg_list = contracted_graph(type_2)
            for c_j in contracted_link_synthesis(type_2):
                answer = conventional_graph(cg_list, c_j)
                answers.extend(answer)
                answer = conventional_graph(cg_list, c_j, 2)
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
        self.assertEqual(8, len(vpoints))
        exprs = vpoints_configure(vpoints, inputs)
        self.assertEqual([
            ('PLAP', 'P0', 'L0', 'a0', 'P1'),
            ('PLLP', 'P2', 'L1', 'L2', 'P1', 'P3'),
            ('PLLP', 'P2', 'L3', 'L4', 'P3', 'P4'),
            ('PLLP', 'P1', 'L5', 'L6', 'P2', 'P5'),
            ('PLLP', 'P5', 'L7', 'L8', 'P4', 'P6'),
            ('PLLP', 'P5', 'L9', 'L10', 'P6', 'P7')
        ], exprs.as_list())
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
            for factors in link_synthesis(nl, nj):
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
            'beta_min': 0.2,
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