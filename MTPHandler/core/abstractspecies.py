import sdRDM

from typing import Dict, Optional, Union
from pydantic import PrivateAttr, field_validator, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .vessel import Vessel


@forge_signature
class AbstractSpecies(sdRDM.DataModel):
    """This object is used to inherit basic attributes common to all species used in the data model."""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: str = element(
        description="None",
        tag="name",
        json_schema_extra=dict(),
    )

    vessel_id: Union[Vessel, str] = element(
        description="None",
        tag="vessel_id",
        json_schema_extra=dict(reference="Vessel.id"),
    )

    init_conc: Optional[float] = element(
        description="None",
        default=None,
        tag="init_conc",
        json_schema_extra=dict(),
    )

    constant: bool = element(
        description="None",
        tag="constant",
        json_schema_extra=dict(),
    )

    unit: Optional[str] = element(
        description="None",
        default=None,
        tag="unit",
        json_schema_extra=dict(),
    )

    uri: Optional[str] = element(
        description="None",
        default=None,
        tag="uri",
        json_schema_extra=dict(),
    )

    creator_id: Optional[str] = element(
        description="None",
        default=None,
        tag="creator_id",
        json_schema_extra=dict(),
    )
    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/MTPHandler"
    )
    _commit: Optional[str] = PrivateAttr(
        default="fce12c40347b8116f04f3d4da2323906c7bf4c7e"
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

    @field_validator("vessel_id")
    def get_vessel_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""
        from .vessel import Vessel

        if isinstance(value, Vessel):
            return value.id
        elif isinstance(value, str):
            return value
        else:
            raise TypeError(
                f"Expected types [Vessel, str] got '{type(value).__name__}' instead."
            )

    @field_validator("vessel_id")
    def get_vessel_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""
        from .vessel import Vessel

        if isinstance(value, Vessel):
            return value.id
        elif isinstance(value, str):
            return value
        else:
            raise TypeError(
                f"Expected types [Vessel, str] got '{type(value).__name__}' instead."
            )
