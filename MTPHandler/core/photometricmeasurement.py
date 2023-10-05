import sdRDM

from typing import List, Optional
from pydantic import Field
from pydantic import Field
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractspecies import AbstractSpecies
from .blankstate import BlankState


@forge_signature
class PhotometricMeasurement(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("photometricmeasurementINDEX"),
        xml="@id",
    )

    wavelength: float = Field(..., description="Wavelength of the measurement")

    wavelength_unit: str = Field(..., description="Unit of the wavelength")

    absorptions: List[float] = Field(
        description="Absorption of the species", multiple=True, default_factory=ListPlus
    )

    blank_states: List[BlankState] = Field(
        description=(
            "List of blank states, referring to the blank state of the species of the"
            " well"
        ),
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_blank_states(
        self,
        species_id: AbstractSpecies,
        was_blanked: bool = False,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'BlankState' to attribute blank_states

        Args:
            id (str): Unique identifier of the 'BlankState' object. Defaults to 'None'.
            species_id (): Reference to species.
            was_blanked (): Whether the species' absorption contribution was subtracted from the absorption signal. Defaults to False
        """
        params = {"species_id": species_id, "was_blanked": was_blanked}
        if id is not None:
            params["id"] = id
        self.blank_states.append(BlankState(**params))
        return self.blank_states[-1]
