import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .well import Well


@forge_signature
class Plate(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("plateINDEX"),
        xml="@id",
    )

    n_rows: Optional[int] = Field(
        default=None,
        description="Number of rows on the plate",
    )

    n_columns: Optional[int] = Field(
        default=None,
        description="Number of columns on the plate",
    )

    temperature: Optional[float] = Field(
        default=None,
        description="Thermostat temperature",
    )

    temperature_unit: Optional[str] = Field(
        default=None,
        description="Unit of the temperature",
    )

    ph: Optional[float] = Field(
        default=None,
        description="pH of the reaction",
    )

    wells: List[Well] = Field(
        description="List of wells on the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_wells(
        self,
        absorption: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        time_unit: Optional[str] = None,
        reaction_volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
        init_conc: Optional[float] = None,
        conc_unit: Optional[str] = None,
        x_position: Optional[int] = None,
        y_position: Optional[int] = None,
        species_id: Optional[str] = None,
        wavelegth: Optional[int] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Well' to attribute wells

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            absorption (): Absorption of the species. Defaults to ListPlus()
            time (): Time of the measurement. Defaults to ListPlus()
            time_unit (): Unit of the time. Defaults to None
            reaction_volume (): Volume of the reaction. Defaults to None
            volume_unit (): Unit of the volume. Defaults to None
            init_conc (): Initial concentration of the species. Defaults to None
            conc_unit (): Concentration unit. Defaults to None
            x_position (): X position of the well on the plate. Defaults to None
            y_position (): Y position of the well on the plate. Defaults to None
            species_id (): ID of the species. Defaults to None
            wavelegth (): Wavelength of the measurement. Defaults to None
        """

        params = {
            "absorption": absorption,
            "time": time,
            "time_unit": time_unit,
            "reaction_volume": reaction_volume,
            "volume_unit": volume_unit,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
            "x_position": x_position,
            "y_position": y_position,
            "species_id": species_id,
            "wavelegth": wavelegth,
        }

        if id is not None:
            params["id"] = id

        self.wells.append(Well(**params))

        return self.wells[-1]
