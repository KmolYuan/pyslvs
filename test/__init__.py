# -*- coding: utf-8 -*-

from unittest import TestCase
from time import time


class TestBase(TestCase):

    def setUp(self) -> None:
        self.t0 = time()

    def tearDown(self) -> None:
        t = time() - self.t0
        print(f"{self.id()} {t * 1000:.04f}ms")
