import sdRDM

from typing import Optional, Union
from pydantic import Field, validator
from sdRDM.base.utils import forge_signature, IDGenerator


from .abstractspecies import AbstractSpecies


@forge_signature
class InitCondition(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("initconditionINDEX"),
        xml="@id",
    )

    species_id: Union[AbstractSpecies, str, None] = Field(
        default=None,
        reference="AbstractSpecies.id",
        description="Reference to species",
    )

    init_conc: Optional[float] = Field(
        default=None,
        description="Initial concentration of the species",
    )

    conc_unit: Optional[str] = Field(
        default=None,
        description="Concentration unit",
    )

    was_blanked: bool = Field(
        description=(
            "Whether the species' absorption contribution was subtracted from the"
            " absorption signal"
        ),
        default=False,
    )

    @validator("species_id")
    def get_species_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""

        from .abstractspecies import AbstractSpecies

        if isinstance(value, AbstractSpecies):
            return value.id
        elif isinstance(value, str):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [AbstractSpecies, str] got '{type(value).__name__}'"
                " instead."
            )
