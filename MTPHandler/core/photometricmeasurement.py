import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from CaliPytion.core import Standard
from .abstractspecies import AbstractSpecies
from .blankstate import BlankState


@forge_signature
class PhotometricMeasurement(sdRDM.DataModel, search_mode="unordered"):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    wavelength: float = element(
        description="Wavelength of the measurement",
        tag="wavelength",
        json_schema_extra=dict(),
    )

    wavelength_unit: str = element(
        description="Unit of the wavelength",
        tag="wavelength_unit",
        json_schema_extra=dict(),
    )

    absorptions: List[float] = element(
        description="Absorption of the species",
        default_factory=ListPlus,
        tag="absorptions",
        json_schema_extra=dict(multiple=True),
    )

    blank_states: List[BlankState] = element(
        description=(
            "List of blank states, referring to the blank state of the species of the"
            " well"
        ),
        default_factory=ListPlus,
        tag="blank_states",
        json_schema_extra=dict(multiple=True),
    )
    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/MTPHandler"
    )
    _commit: Optional[str] = PrivateAttr(
        default="e87642023bceb2ac5538980efc1e78fd8e7164b4"
    )
    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                (isinstance(i, _Element) for i in value)
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)
        return self

    def add_to_blank_states(
        self,
        species_id: AbstractSpecies,
        contributes_to_signal: bool = True,
        id: Optional[str] = None,
    ) -> BlankState:
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

    def to_concentration(
        self, standard: Standard, ignore_blank_status: bool = False, **kwargs
    ) -> list[float]:
        if not standard.model_result.was_fitted:
            raise ValueError(
                f"Standard {standard.id} was not fitted for species"
                f" {standard.species_id}"
            )

        if not self.wavelength == standard.wavelength:
            raise ValueError(
                f"Standard at {standard.wavelength} nm not applicable for well"
                f" {self.id} measured at {self.wavelength}"
            )

        if not ignore_blank_status and not self.is_blanked_for(standard.species_id):
            raise ValueError(
                f"Well {self.id} was not blanked for species {standard.species_id}"
                f"{self.blank_states}"
            )

        return standard.model_result.calculate(self.absorptions, **kwargs)
