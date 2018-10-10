# -*- coding: utf-8 -*-

from typing import Any
import numpy as np


class Verification:

    """Verification function class base.

    The 'verify' module should be loaded when using sub-class.
    """

    def __call__(self, v: np.ndarray) -> float:
        """Calculate the fitness.

        Usage:
        f = MyVerification()
        fitness = f(chromosome.v)
        """
        ...

    def result(self, v: np.ndarray) -> Any:
        """Show the result."""
        ...
