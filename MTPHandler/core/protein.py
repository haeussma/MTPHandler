from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.datatypes import Identifier
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict


class Protein(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Description of a protein species that might be present in the wells of the plate."""

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

    sequence: Optional[str] = element(
        description="Amino acid sequence of the protein",
        default=None,
        tag="sequence",
        json_schema_extra=dict(),
    )

    organism: Optional[str] = element(
        description="Organism the protein originates from",
        default=None,
        tag="organism",
        json_schema_extra=dict(),
    )

    organism_tax_id: Optional[Identifier] = element(
        description="NCBI taxonomy ID of the organism",
        default=None,
        tag="organism_tax_id",
        json_schema_extra=dict(),
    )

    references: List[Identifier] = element(
        description="List of references to the protein",
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
        default="e334b0b111f8283b76d4ff24a987827a4cff7116"
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
