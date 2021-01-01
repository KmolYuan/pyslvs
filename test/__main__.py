# -*- coding: utf-8 -*-

"""This module will test the functions of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import defaultTestLoader, TextTestRunner, TestSuite
from test_metaheuristics import AlgorithmTest

if __name__ == '__main__':
    TextTestRunner().run(TestSuite([
        defaultTestLoader.discover('test'),
        defaultTestLoader.loadTestsFromTestCase(AlgorithmTest),
    ]))
