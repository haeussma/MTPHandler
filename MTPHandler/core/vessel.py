import sdRDM

from typing import Dict, Optional
from pydantic import PositiveFloat, PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class Vessel(sdRDM.DataModel, search_mode="unordered"):
    """This object describes vessels in which the experiment has been carried out. These can include any type of vessel used in biocatalytic experiments."""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    name: str = element(
        description="Name of the used vessel.",
        tag="name",
        json_schema_extra=dict(template_alias="Name"),
    )

    volume: PositiveFloat = element(
        description="Volumetric value of the vessel.",
        tag="volume",
        json_schema_extra=dict(template_alias="Volume value"),
    )

    unit: str = element(
        description="Volumetric unit of the vessel.",
        tag="unit",
        json_schema_extra=dict(template_alias="Volume unit"),
    )

    constant: bool = element(
        description="Whether the volume of the vessel is constant or not.",
        default=True,
        tag="constant",
        json_schema_extra=dict(),
    )

    uri: Optional[str] = element(
        description="URI of the vessel.",
        default=None,
        tag="uri",
        json_schema_extra=dict(),
    )

    creator_id: Optional[str] = element(
        description="Unique identifier of the author.",
        default=None,
        tag="creator_id",
        json_schema_extra=dict(),
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
