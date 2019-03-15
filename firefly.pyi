# -*- coding: utf-8 -*-

from typing import (
    Tuple,
    List,
    Dict,
    Callable,
    Optional,
    Any,
)
from .verify import Verification


class Firefly:

    """Algorithm class."""

    def __init__(
        self,
        func: Verification,
        settings: Dict[str, Any],
        progress_fun: Optional[Callable[[int, str], None]] = None,
        interrupt_fun: Optional[Callable[[], bool]] = None
    ):
        """
        settings = {
            'n',
            'alpha',
            'beta_min',
            'beta0',
            'gamma',
            'max_gen', 'min_fit' or 'max_time',
            'report'
        }
        """
        ...

    def run(self) -> Tuple[Any, List[Tuple[int, float, float]]]:
        ...
