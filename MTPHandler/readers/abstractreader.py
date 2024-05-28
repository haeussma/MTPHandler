from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate


class AbstractReader(ABC):
    """
    AbstractReader is an abstract base class that defines the
    interface for reading and processing data from a file containing
    measurement data from a microtiter plate.

    Attributes:
        path (str): The path to the measurement file.
    """

    def __init__(self, path: str):
        self.path = path

    @abstractmethod
    def read(self) -> Plate:
        pass
