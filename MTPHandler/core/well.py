from typing import Dict, List, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.datatypes import Unit
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from .blankstate import BlankState
from .initcondition import InitCondition
from .photometricmeasurement import PhotometricMeasurement


class Well(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Description of a well on the plate."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    x_pos: int = element(
        description="X position of the well on the plate",
        tag="x_pos",
        json_schema_extra=dict(),
    )

    y_pos: int = element(
        description="Y position of the well on the plate",
        tag="y_pos",
        json_schema_extra=dict(),
    )

    ph: Optional[float] = element(
        description="pH of the reaction",
        default=None,
        tag="ph",
        json_schema_extra=dict(),
    )

    init_conditions: List[InitCondition] = element(
        description="List of initial conditions of different species",
        default_factory=ListPlus,
        tag="init_conditions",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    measurements: List[PhotometricMeasurement] = element(
        description="List of photometric measurements",
        default_factory=ListPlus,
        tag="measurements",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    volume: Optional[float] = element(
        description="Volume of the reaction",
        default=None,
        tag="volume",
        json_schema_extra=dict(),
    )

    volume_unit: Optional[Unit] = element(
        description="Unit of the volume",
        default=None,
        tag="volume_unit",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/MTPHandler"
    )
    _commit: Optional[str] = PrivateAttr(
        default="c1ace6fef38751e79952a37557221f0540ca9d77"
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

    def add_to_init_conditions(
        self,
        species_id: str,
        init_conc: float,
        conc_unit: Unit,
        id: Optional[str] = None,
        **kwargs,
    ) -> InitCondition:
        """
        This method adds an object of type 'InitCondition' to attribute init_conditions

        Args:
            id (str): Unique identifier of the 'InitCondition' object. Defaults to 'None'.
            species_id (): Reference to species.
            init_conc (): Initial concentration of the species.
            conc_unit (): Concentration unit.
        """

        params = {
            "species_id": species_id,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
        }

        if id is not None:
            params["id"] = id

        obj = InitCondition(**params)

        self.init_conditions.append(obj)

        return self.init_conditions[-1]

    def add_to_measurements(
        self,
        wavelength: float,
        wavelength_unit: Unit,
        absorption: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        time_unit: Optional[Unit] = None,
        blank_states: List[BlankState] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> PhotometricMeasurement:
        """
        This method adds an object of type 'PhotometricMeasurement' to attribute measurements

        Args:
            id (str): Unique identifier of the 'PhotometricMeasurement' object. Defaults to 'None'.
            wavelength (): Wavelength of the measurement.
            wavelength_unit (): Unit of the wavelength.
            absorption (): Absorption of the species. Defaults to ListPlus()
            time (): Time of the measurement. Defaults to ListPlus()
            time_unit (): Unit of the time. Defaults to None
            blank_states (): List of blank states, referring to the blank state of the species of the well. Defaults to ListPlus()
        """

        params = {
            "wavelength": wavelength,
            "wavelength_unit": wavelength_unit,
            "absorption": absorption,
            "time": time,
            "time_unit": time_unit,
            "blank_states": blank_states,
        }

        if id is not None:
            params["id"] = id

        obj = PhotometricMeasurement(**params)

        self.measurements.append(obj)

        return self.measurements[-1]

    def _contains_species(self, species_id: str) -> bool:
        for condition in self.init_conditions:
            if condition.species_id == species_id:
                return True

        return False

    def _get_species_condition(self, species_id: str) -> InitCondition:
        for condition in self.init_conditions:
            if condition.species_id == species_id:
                return condition

        raise ValueError(f"Species {species_id} not found in well {self.id}")

    def get_measurement(self, wavelength: float) -> PhotometricMeasurement:
        for measurement in self.measurements:
            if measurement.wavelength == wavelength:
                return measurement

        raise ValueError(f"No measurement at {wavelength} nm found for well {self.id}")
