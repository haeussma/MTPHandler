import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractspecies import AbstractSpecies
from .blankstate import BlankState
from .initcondition import InitCondition
from .photometricmeasurement import PhotometricMeasurement


@forge_signature
class Well(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("wellINDEX"),
        xml="@id",
    )

    ph: float = Field(
        ...,
        description="pH of the reaction",
    )

    x_position: int = Field(
        ...,
        description="X position of the well on the plate",
    )

    y_position: int = Field(
        ...,
        description="Y position of the well on the plate",
    )

    init_conditions: List[InitCondition] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="List of initial conditions of different species",
    )

    measurements: List[PhotometricMeasurement] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="List of photometric measurements",
    )

    volume: Optional[float] = Field(
        default=None,
        description="Volume of the reaction",
    )

    volume_unit: Optional[str] = Field(
        default=None,
        description="Unit of the volume",
    )

    def add_to_init_conditions(
        self,
        species_id: AbstractSpecies,
        init_conc: float,
        conc_unit: str,
        id: Optional[str] = None,
    ) -> None:
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
    ) -> None:
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
