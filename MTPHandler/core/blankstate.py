import sdRDM

from typing import Dict, Optional, Union
from pydantic import PrivateAttr, field_validator, model_validator, validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .abstractspecies import AbstractSpecies


@forge_signature
class BlankState(sdRDM.DataModel):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    species_id: Union[AbstractSpecies, str] = element(
        description="Reference to species",
        tag="species_id",
        json_schema_extra=dict(reference="AbstractSpecies.id"),
    )

    contributes_to_signal: bool = element(
        description=(
            "Whether the species' absorption contributes to the absorption signal"
        ),
        default=True,
        tag="contributes_to_signal",
        json_schema_extra=dict(),
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

    @field_validator("species_id")
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
