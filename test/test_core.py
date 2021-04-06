# -*- coding: utf-8 -*-

"""Pyslvs core module test."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from math import sqrt, radians
from pyslvs import (
    Coord, SolverSystem, pxy, ppp, plap, pllp, plpp, palp, expr_solving,
    t_config, parse_vpoints, example_list,
)
from . import TestBase


class CoreTest(TestBase):

    def test_pxy(self):
        """Test for pxy function."""
        coord = pxy(Coord(80, 90), 40, -20)
        self.assertAlmostEqual(120, coord.x)
        self.assertAlmostEqual(70, coord.y)

    def test_ppp(self):
        """Test for ppp function."""
        coord = ppp(Coord(0, 0), Coord(0, 90), Coord(90, 0))
        self.assertAlmostEqual(90, coord.x)
        self.assertAlmostEqual(90, coord.y)

    def test_plap(self):
        """Test for plap function."""
        coord = plap(Coord(0, 0), 50 * sqrt(2), radians(45), Coord(50, 0))
        self.assertAlmostEqual(50, coord.x)
        self.assertAlmostEqual(50, coord.y)

    def test_pllp(self):
        """Test for pllp function."""
        c1 = Coord(-30, 0)
        c2 = Coord(30, 0)
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
        coord = plpp(Coord(0, 0), sqrt(5), Coord(0, -3), Coord(3 / 2, 0))
        self.assertAlmostEqual(2, coord.x)
        self.assertAlmostEqual(1, coord.y)

    def test_palp(self):
        """Test for palp function."""
        coord = palp(Coord(0, 0), radians(15), 20, Coord(60, 10))
        self.assertAlmostEqual(42.253221, coord.x, 6)
        self.assertAlmostEqual(19.222356, coord.y, 6)

    def test_solving(self):
        """Test triangular formula solving.

        + Test for PMKS parser.
        + Test data collecting function.
        + Test expression solving function.
        """

        def test_case(name: str):
            expr, inputs = example_list(name)
            vpoints = parse_vpoints(expr)
            exprs = t_config(vpoints, inputs)
            result = expr_solving(exprs, vpoints, {pair: 0. for pair in inputs})
            return result[-1]

        x, y = test_case("Jansen's linkage (Single)")
        self.assertAlmostEqual(-43.170055, x, 6)
        self.assertAlmostEqual(-91.753226, y, 6)
        x, y = test_case("Crank slider (RP joint)")
        self.assertAlmostEqual(103.801126, x, 6)
        self.assertAlmostEqual(78.393173, y, 6)
        x, y = test_case("Parallel Linkage")
        self.assertAlmostEqual(200, x, 6)
        self.assertAlmostEqual(0, y, 6)
        # TODO: New test case for Inverted slider

    def test_solving_bfgs(self):
        """Test Sketch Solve kernel."""
        expr, _ = example_list("Jansen's linkage (Single)")
        system = SolverSystem(parse_vpoints(expr), {(0, 1): 0.})
        result = system.solve()
        x, y = result[7]
        self.assertAlmostEqual(-43.170055, x, 6)
        self.assertAlmostEqual(-91.753226, y, 6)
        # Test if angle value changed
        system.set_inputs({(0, 1): 45.})
        result = system.solve()
        x, y = result[7]
        self.assertAlmostEqual(-24.406394, x, 6)
        self.assertAlmostEqual(-91.789596, y, 6)
        # Test if link length changed
        system.set_data({(0, 1): 16.})
        result = system.solve()
        x, y = result[7]
        self.assertAlmostEqual(-24.117994, x, 6)
        self.assertAlmostEqual(-91.198072, y, 6)
