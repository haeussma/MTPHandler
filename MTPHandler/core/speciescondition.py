from re import T
import sdRDM

from typing import Optional, Union
from pydantic import Field, validator
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species
from .speciestype import SpeciesType


@forge_signature
class SpeciesCondition(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("speciesconditionINDEX"),
        xml="@id",
    )

    species_type: Union[Species, SpeciesType, None] = Field(
        default=None,
        reference="Species.type",
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

    @validator("species_type", pre=True)
    def get_species_type_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""

        from .species import Species

        if isinstance(value, Species):
            return value.type
        elif isinstance(value, SpeciesType):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [Species, SpeciesType] got '{type(value).__name__}'"
                " instead."
            )
