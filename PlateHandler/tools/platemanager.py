import sdRDM

from typing import Optional, List
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator

from PlateHandler.core import Plate
from PlateHandler.tools.spectramax_reader import read_spectramax


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

    def assign_concentration(
            self,
            id: str,
            concentration: float,
            concentration_unit: str,
    ):
        wells = self.plate.get("wells", "id", id.upper())[0]

        for well in wells:
            well.init_conc = concentration
            well.conc_unit = concentration_unit

    def assign_concentration_series(
            self,
            concentration: List[float],
            concentration_unit: str,
            to_rows: List[str] = None,
            to_columns: List[int] = None,
            wavelength: int = None
    ):
        if (to_rows and not to_columns) or (not to_rows and to_columns):
            raise AttributeError(
                "Specify either 'to_rows' or 'to_columns'."
            )

        if to_rows:
            pass
