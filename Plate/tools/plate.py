import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from astropy.units import UnitBase

from ..core.well import Well
from ..core.abstractspecies import AbstractSpecies


@forge_signature
class Plate(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("plateINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the plate",
    )

    n_rows: Optional[int] = Field(
        default=None,
        description="Number of rows on the plate",
    )

    n_columns: Optional[int] = Field(
        default=None,
        description="Number of columns on the plate",
    )

    wells: List[Well] = Field(
        description="List of wells on the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_wells(
        self,
        id: str,
        reaction_volume: Optional[float] = None,
        unit: Optional[UnitBase] = None,
        x_position: Optional[int] = None,
        y_position: Optional[int] = None,
        species: List[AbstractSpecies] = ListPlus(),
    ) -> None:
        """
        This method adds an object of type 'Well' to attribute wells

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            reaction_volume (): Volume of the reaction. Defaults to None
            unit (): Unit of the volume. Defaults to None
            x_position (): X position of the well on the plate. Defaults to None
            y_position (): Y position of the well on the plate. Defaults to None
            species (): List of species present in the well. Defaults to ListPlus()
        """

        params = {
            "reaction_volume": reaction_volume,
            "unit": unit,
            "x_position": x_position,
            "y_position": y_position,
            "species": species,
        }

        if id is not None:
            params["id"] = id

        self.wells.append(Well(**params))

        return self.wells[-1]
