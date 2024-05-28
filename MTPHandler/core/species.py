from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.datatypes import Identifier
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class Species(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Description of a species that might be present in the wells of the plate."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    name: Optional[str] = element(
        description="Name of the species",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    references: List[Identifier] = element(
        description="List of references to the species",
        default_factory=ListPlus,
        tag="references",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/MTPHandler"
    )
    _commit: Optional[str] = PrivateAttr(
        default="7c8baa9049c8c8a4520b326b921931efc637101b"
    )

    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                isinstance(i, _Element) for i in value
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)

        return self
