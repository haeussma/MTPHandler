import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Well(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("wellINDEX"),
        xml="@id",
    )

    absorption: List[float] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Absorption of the species",
    )

    time: List[float] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Time of the measurement",
    )

    time_unit: Optional[str] = Field(
        default=None,
        description="Unit of the time",
    )

    reaction_volume: Optional[float] = Field(
        default=None,
        description="Volume of the reaction",
    )

    volume_unit: Optional[str] = Field(
        default=None,
        description="Unit of the volume",
    )

    init_conc: Optional[float] = Field(
        default=None,
        description="Initial concentration of the species",
    )

    conc_unit: Optional[str] = Field(
        default=None,
        description="Concentration unit",
    )

    x_position: Optional[int] = Field(
        default=None,
        description="X position of the well on the plate",
    )

    y_position: Optional[int] = Field(
        default=None,
        description="Y position of the well on the plate",
    )

    species_id: Optional[str] = Field(
        default=None,
        description="ID of the species",
    )

    wavelegth: Optional[int] = Field(
        default=None,
        description="Wavelength of the measurement",
    )
