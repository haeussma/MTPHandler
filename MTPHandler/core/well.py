import sdRDM

from typing import Dict, List, Optional
from pydantic import PrivateAttr, model_validator
from uuid import uuid4
from pydantic_xml import attr, element
from lxml.etree import _Element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict
from .blankstate import BlankState
from .initcondition import InitCondition
from .photometricmeasurement import PhotometricMeasurement
from .abstractspecies import AbstractSpecies


@forge_signature
class Well(sdRDM.DataModel):
    """"""

    id: Optional[str] = attr(
        name="id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
        xml="@id",
    )

    ph: float = element(
        description="pH of the reaction",
        tag="ph",
        json_schema_extra=dict(),
    )

    x_position: int = element(
        description="X position of the well on the plate",
        tag="x_position",
        json_schema_extra=dict(),
    )

    y_position: int = element(
        description="Y position of the well on the plate",
        tag="y_position",
        json_schema_extra=dict(),
    )

    init_conditions: List[InitCondition] = element(
        description="List of initial conditions of different species",
        default_factory=ListPlus,
        tag="init_conditions",
        json_schema_extra=dict(multiple=True),
    )

    measurements: List[PhotometricMeasurement] = element(
        description="List of photometric measurements",
        default_factory=ListPlus,
        tag="measurements",
        json_schema_extra=dict(multiple=True),
    )

    volume: Optional[float] = element(
        description="Volume of the reaction",
        default=None,
        tag="volume",
        json_schema_extra=dict(),
    )

    volume_unit: Optional[str] = element(
        description="Unit of the volume",
        default=None,
        tag="volume_unit",
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

    def add_to_init_conditions(
        self,
        species_id: AbstractSpecies,
        init_conc: float,
        conc_unit: str,
        id: Optional[str] = None,
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
        self.init_conditions.append(InitCondition(**params))
        return self.init_conditions[-1]

    def add_to_measurements(
        self,
        wavelength: float,
        wavelength_unit: str,
        absorptions: List[float] = ListPlus(),
        blank_states: List[BlankState] = ListPlus(),
        id: Optional[str] = None,
    ) -> PhotometricMeasurement:
        """
        This method adds an object of type 'PhotometricMeasurement' to attribute measurements

        Args:
            id (str): Unique identifier of the 'PhotometricMeasurement' object. Defaults to 'None'.
            wavelength (): Wavelength of the measurement.
            wavelength_unit (): Unit of the wavelength.
            absorptions (): Absorption of the species. Defaults to ListPlus()
            blank_states (): List of blank states, referring to the blank state of the species of the well. Defaults to ListPlus()
        """
        params = {
            "wavelength": wavelength,
            "wavelength_unit": wavelength_unit,
            "absorptions": absorptions,
            "blank_states": blank_states,
        }
        if id is not None:
            params["id"] = id
        self.measurements.append(PhotometricMeasurement(**params))
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
