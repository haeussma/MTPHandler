import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


from .speciestype import SpeciesType


@forge_signature
class Species(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("speciesINDEX"),
        xml="@id",
    )

    type: SpeciesType = Field(
        ...,
        description="Type of the species",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the species",
    )
