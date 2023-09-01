import sdRDM

from typing import Optional, List
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator

from Plate.core import Plate
from Plate.tools.spectramax_reader import read_spectramax


@forge_signature
class PlateManager(sdRDM.DataModel):

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("platemanagerINDEX"),
        xml="@id",
    )

    plate: Plate = Field(
        default=None,
        description="Plate data",
    )

    @classmethod
    def from_file(
        cls,
        path: str,
        time: List[float],
        time_unit: str,
        ph: float = None,
        temperature: float = None,
        temperature_unit: str = None,
    ):
        return cls(plate=read_spectramax(path, time, time_unit, ph, temperature, temperature_unit))
