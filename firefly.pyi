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
            'betaMin',
            'beta0',
            'gamma',
            'maxGen', 'minFit' or 'maxTime',
            'report'
        }
        """
        ...

    def run(self) -> Tuple[Any, List[Tuple[int, float, float]]]:
        ...
