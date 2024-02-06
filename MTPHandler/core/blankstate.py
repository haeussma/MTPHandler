import sdRDM

from typing import Optional, Union
from pydantic import Field, validator
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractspecies import AbstractSpecies


@forge_signature
class BlankState(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("blankstateINDEX"),
        xml="@id",
    )

    species_id: Union[AbstractSpecies, str] = Field(
        ...,
        reference="AbstractSpecies.id",
        description="Reference to species",
    )

    contributes_to_signal: bool = Field(
        description=(
            "Whether the species' absorption contributes to the absorption signal"
        ),
        default=True,
    )

    @validator("species_id")
    def get_species_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""
        from .abstractspecies import AbstractSpecies

        if isinstance(value, AbstractSpecies):
            return value.id
        elif isinstance(value, str):
            return value
        else:
            raise TypeError(
                f"Expected types [AbstractSpecies, str] got '{type(value).__name__}'"
                " instead."
            )
