# -*- coding: utf-8 -*-

"""Pyslvs graph module test."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from pyslvs.graph import (
    link_assortment,
    contracted_link_assortment,
    Graph,
    is_planar,
    external_loop_layout,
    conventional_graph,
    contracted_graph,
    link_synthesis,
    contracted_link_synthesis,
)
from . import TestBase

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


class GraphTest(TestBase):

    def test_graph_basic(self):
        """Test 'graph' libraries."""
        g1 = Graph([(0, 1), (0, 4), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)])
        self.assertFalse(g1.is_degenerate())
        self.assertTrue(g1.is_connected())
        self.assertEqual(1, g1.dof())

    def test_graph_isomorphic(self):
        g1 = Graph([(0, 1), (0, 4), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)])
        g2 = Graph([(0, 2), (0, 4), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)])
        g3 = Graph([(0, 1), (0, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5)])
        self.assertTrue(g1.is_isomorphic(g2))
        self.assertFalse(g1.is_isomorphic(g3))

    def test_graph_planar(self):
        g1 = Graph([
            (0, 1), (2, 7), (1, 5), (1, 6), (3, 6), (0, 4), (3, 7), (2, 5),
            (3, 4), (0, 2),
        ])
        self.assertTrue(is_planar(g1))
        self.assertEqual([4, 4], link_assortment(g1))
        self.assertEqual([4, 0, 0, 0], contracted_link_assortment(g1))

    def test_graph_degenerate(self):
        g1 = Graph([(0, 1), (0, 2), (2, 1), (0, 3), (3, 4), (4, 5), (5, 1)])
        self.assertTrue(g1.is_degenerate())

    def test_graph_duplicate(self):
        g1 = Graph([(0, 1), (1, 2), (2, 3), (0, 3)])
        g2 = g1.duplicate([2, 3], 1)
        self.assertEqual(set(g2.edges), {
            (0, 1), (1, 2), (4, 5), (1, 4), (2, 3), (0, 5), (0, 3),
        })

    def test_graph_loop(self):
        g1 = Graph([
            (0, 1), (4, 10), (1, 3), (2, 9), (5, 6), (4, 5), (5, 7), (8, 10),
            (1, 8), (9, 11), (3, 6), (0, 4), (3, 7), (2, 5), (0, 2), (4, 11),
        ])
        pos = external_loop_layout(g1, True)
        self.assertEqual(set(g1.vertices), set(pos))

    def test_graph_degree_code(self):
        g1 = Graph([
            (0, 1), (0, 2), (0, 3), (0, 5), (2, 3), (3, 4), (2, 5), (3, 5),
            (4, 6), (3, 6),
        ])
        g2 = Graph([
            (0, 1), (0, 4), (0, 2), (0, 6), (0, 3), (1, 4), (2, 6), (2, 5),
            (3, 6), (2, 3),
        ])
        self.assertTrue(g1.is_isomorphic(g2))
        self.assertEqual(2057732, g1.degree_code())
        self.assertEqual(g1.degree_code(), g2.degree_code())
        for code, edges in DEGREE_CODE_TABLE:
            self.assertEqual(code, Graph(edges).degree_code())

    def test_atlas(self):
        """Test 'atlas' libraries."""
        answers = []
        # Test assortment [4]
        type_0 = [4]
        cg_list = contracted_graph(type_0)
        for c_j in contracted_link_synthesis(type_0):
            answer = conventional_graph(cg_list, c_j)
            answers.extend(answer)
        self.assertEqual(1, len(answers))
        answers.clear()
        # Test assortment [4, 2]
        type_1 = [4, 2]
        cg_list = contracted_graph(type_1)
        for c_j in contracted_link_synthesis(type_1):
            answer = conventional_graph(cg_list, c_j)
            answers.extend(answer)
        self.assertEqual(2, len(answers))
        answers.clear()
        # Test assortment [4, 4, 0], [5, 2, 1], [6, 0, 2]
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

    def test_number_synthesis(self):
        """Test Number Synthesis function."""
        for nl, nj in [(4, 4), (6, 7), (8, 9), (10, 12)]:
            for factors in link_synthesis(nl, nj):
                count = 0
                for i, factor in enumerate(factors):
                    count += factor * (i + 2)
                self.assertEqual(nj, int(count / 2))
