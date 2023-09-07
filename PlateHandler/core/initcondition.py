import sdRDM

from typing import Optional, Union
from pydantic import Field, validator
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species


@forge_signature
class InitCondition(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("initconditionINDEX"),
        xml="@id",
    )

    species: Union[Species, str, None] = Field(
        default=None,
        reference="Species.id",
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

    @validator("species")
    def get_species_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""

        from .species import Species

        if isinstance(value, Species):
            return value.id
        elif isinstance(value, str):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [Species, str] got '{type(value).__name__}' instead."
            )
