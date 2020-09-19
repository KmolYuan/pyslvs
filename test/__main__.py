# -*- coding: utf-8 -*-

"""This module will test the functions of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from unittest import defaultTestLoader, TextTestRunner

if __name__ == '__main__':
    TextTestRunner().run(defaultTestLoader.discover('test'))
