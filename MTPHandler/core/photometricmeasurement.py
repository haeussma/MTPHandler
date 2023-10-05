import sdRDM

from typing import List, Optional
from pydantic import Field
from pydantic import Field
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .blankstate import BlankState
from .abstractspecies import AbstractSpecies


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
        contributes_to_signal: bool = True,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'BlankState' to attribute blank_states

        Args:
            id (str): Unique identifier of the 'BlankState' object. Defaults to 'None'.
            species_id (): Reference to species.
            contributes_to_signal (): Whether the species' absorption contributes to the absorption signal. Defaults to True
        """
        params = {
            "species_id": species_id,
            "contributes_to_signal": contributes_to_signal,
        }
        if id is not None:
            params["id"] = id
        self.blank_states.append(BlankState(**params))
        return self.blank_states[-1]

    def is_blanked_for(self, species_id: str) -> bool:
        """Checks if the measurement is blanked for a given species."""

        if species_id not in [state.species_id for state in self.blank_states]:
            raise ValueError(f"Species {species_id} is not present in this well.")

        target_contributes = [
            state.contributes_to_signal
            for state in self.blank_states
            if state.species_id == species_id
        ][0]

        others_contribute = [
            state.contributes_to_signal
            for state in self.blank_states
            if state.species_id != species_id
        ]

        if target_contributes and not any(others_contribute):
            return True

        return False

    def species_contibutes(self, species_id: str) -> bool:
        species_contributes = [
            state.contributes_to_signal
            for state in self.blank_states
            if state.species_id == species_id
        ][0]

        return species_contributes

    def get_blank_state(self, species_id: str) -> BlankState:
        for state in self.blank_states:
            if state.species_id == species_id:
                return state

        raise ValueError(f"Species {species_id} is not present in this well.")
