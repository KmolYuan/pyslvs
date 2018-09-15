# -*- coding: utf-8 -*-

"""This module will test the functions of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2018"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

import unittest
from unittest import TestCase
from typing import Tuple, List

# For necessary testing modules.
from math import sqrt, radians, isclose
from pmks import VPoint
from bfgs import vpoint_solving
from tinycadlib import (
    Coordinate,
    PLAP,
    PLLP,
    PLPP,
    PXY,
    expr_solving,
    data_collecting,
)
from planarlinkage import Planar
from rga import Genetic
from firefly import Firefly
from de import DiffertialEvolution
from number import number_synthesis
from atlas import topo, Graph
from triangulation import vpoints_configure
from _parser import parse_vpoints
from examples import example_list


class CoreTest(TestCase):
    
    """Testing Cython libs."""
    
    def vpoints_object(self) -> Tuple[List[VPoint], Tuple[Tuple[int, int]]]:
        """Example: Jansen's linkage (Single)."""
        expr, inputs = example_list["Jansen's linkage (Single)"]
        return parse_vpoints(expr), inputs
    
    def planar_object(self) -> Planar:
        """Test-used mechanism for algorithm."""
        return Planar({
            'Driver': {'P0': (-70, -70, 50)},
            'Follower': {'P1': (70, -70, 50)},
            'Target': {'P4': [
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
            ]},
            'Expression':
                "PLAP[P0,L0,a0](P2);"
                "PLLP[P2,L1,L2,P1](P3);"
                "PLLP[P2,L3,L4,P3](P4)",
            'constraints': [('P0', 'P1', 'P2', 'P3')],
            'upper': [100., 100., 100., 100., 100., 360.],
            'lower': [5., 5., 5., 5., 5., 0.],
        })
    
    def test_plap(self):
        """Test for PLAP function."""
        A = Coordinate(0, 0)
        B = Coordinate(50, 0)
        x, y = PLAP(A, 50 * sqrt(2), radians(45), B)
        self.assertTrue(isclose(x, 50))
        self.assertTrue(isclose(y, 50))
    
    def test_pllp(self):
        """Test for PLLP function."""
        A = Coordinate(-30, 0)
        B = Coordinate(30, 0)
        x, y = PLLP(A, 50, 50, B)
        self.assertTrue(isclose(x, 0))
        self.assertTrue(isclose(y, 40))
        x, y = PLLP(A, 30, 30, B)
        self.assertTrue(isclose(x, 0))
        self.assertTrue(isclose(y, 0))
        x, y = PLLP(A, 90, 30, B)
        self.assertTrue(isclose(x, 60))
        self.assertTrue(isclose(y, 0))
    
    def test_plpp(self):
        """Test for PLPP function."""
        A = Coordinate(0, 0)
        B = Coordinate(0, -3)
        C = Coordinate(3/2, 0)
        x, y = PLPP(A, sqrt(5), B, C)
        self.assertTrue(isclose(x, 2))
        self.assertTrue(isclose(y, 1))
    
    def test_pxy(self):
        A = Coordinate(80, 90)
        x, y = PXY(A, 40, -20)
        self.assertTrue(isclose(x, 120))
        self.assertTrue(isclose(y, 70))
    
    def test_topologic(self):
        """Testing 'topologic' libraries.
        
        + 'topo' function.
        + 'Graph' class.
        """
        G = Graph([(0, 1), (0, 4), (1, 5), (2, 3), (2, 4), (3, 5), (4, 5)])
        H = Graph([(0, 2), (0, 4), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)])
        I = Graph([(0, 1), (0, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5)])
        self.assertTrue(G.is_isomorphic(H))
        self.assertFalse(G.is_isomorphic(I))
        answer, time = topo([4, 2], degenerate=True)
        self.assertEqual(len(answer), 2)
    
    def test_solving(self):
        """Test triangular formula solving.
        
        + Test for PMKS parser.
        + Test data collecting function.
        + Test expression solving function.
        """
        vpoints, inputs = self.vpoints_object()
        self.assertTrue(len(vpoints) == 8)
        exprs = vpoints_configure(vpoints, inputs)
        mapping = {n: f'P{n}' for n in range(len(vpoints))}
        data_dict, dof = data_collecting(exprs, mapping, vpoints)
        for link, link_length in (
            ('L0', 15.002083188677497),
            ('L1', 41.50187586121861),
            ('L2', 49.9949057404852),
            ('L3', 40.09650982317538),
            ('L4', 55.80253220060896),
            ('L5', 61.90525179659639),
            ('L6', 39.302800154696364),
            ('L7', 36.69767567571548),
            ('L8', 39.395233214184685),
            ('L9', 48.995886562037015),
            ('L10', 65.69940106271898),
        ):
            self.assertTrue(isclose(data_dict[link], link_length))
        self.assertEqual(dof, 1)
        result = expr_solving(exprs, mapping, vpoints, [0.])
        x, y = result[-1]
        self.assertTrue(isclose(x, -43.17005515543241))
        self.assertTrue(isclose(y, -91.75322590542523))
    
    def test_bfgs(self):
        """Test Sketch Solve kernel."""
        vpoints, inputs = self.vpoints_object()
        result = vpoint_solving(vpoints, [(0, 1, 0.)])
        x, y = result[-1]
        self.assertTrue(isclose(round(x, 2), -43.17))
        self.assertTrue(isclose(round(y, 2), -91.75))
    
    def test_number_synthesis(self):
        """Test Number Synthesis function."""
        for NL, NJ in [(4, 4), (6, 7), (8, 9), (10, 12)]:
            for factors in number_synthesis(NL, NJ):
                count = 0
                for i, factor in enumerate(factors):
                    count += factor * (i + 2)
                self.assertEqual(int(count / 2), NJ)
    
    def test_algorithm_rga(self):
        """Real-coded genetic algorithm."""
        fun1 = Genetic(self.planar_object(), {
            'maxTime': 1,
            'report': 10,
            # Genetic
            'nPop': 500,
            'pCross': 0.95,
            'pMute': 0.05,
            'pWin': 0.95,
            'bDelta': 5.,
        })
        fun1.run()
    
    def test_algorithm_firefly(self):
        """Firefly algorithm."""
        fun2 = Firefly(self.planar_object(), {
            'maxTime': 1,
            'report': 10,
            # Firefly
            'n': 80,
            'alpha': 0.01,
            'betaMin': 0.2,
            'gamma': 1.,
            'beta0': 1.,
        })
        fun2.run()
    
    def test_algorithm_de(self):
        """Differential evolution."""
        fun3 = DiffertialEvolution(self.planar_object(), {
            'maxTime': 1,
            'report': 10,
            # DE
            'strategy': 1,
            'NP': 400,
            'F': 0.6,
            'CR': 0.9,
        })
        fun3.run()


if __name__ == '__main__':
    unittest.main()
