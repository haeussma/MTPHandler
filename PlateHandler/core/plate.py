import sdRDM

from typing import Optional, Union, List
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .initcondition import InitCondition
from .well import Well
from .abstractspecies import AbstractSpecies
from .reactant import Reactant
from .protein import Protein


@forge_signature
class Plate(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("plateINDEX"),
        xml="@id",
    )

    n_rows: Optional[int] = Field(
        default=None,
        description="Number of rows on the plate",
    )

    n_columns: Optional[int] = Field(
        default=None,
        description="Number of columns on the plate",
    )

    temperature: Optional[float] = Field(
        default=None,
        description="Thermostat temperature",
    )

    temperature_unit: Optional[str] = Field(
        default=None,
        description="Unit of the temperature",
    )

    ph: Optional[float] = Field(
        default=None,
        description="pH of the reaction",
    )

    wells: List[Well] = Field(
        description="List of wells on the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    measured_wavelengths: List[int] = Field(
        description="Measured wavelengths in nm",
        default_factory=ListPlus,
        multiple=True,
    )

    species: Union[AbstractSpecies, Protein, Reactant, None] = Field(
        default=None,
        description="List of species present in wells of the plate",
    )

    def add_to_wells(
        self,
        absorption: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        time_unit: Optional[str] = None,
        reaction_volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
        init_conditions: List[InitCondition] = ListPlus(),
        x_position: Optional[int] = None,
        y_position: Optional[int] = None,
        wavelength: Optional[int] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Well' to attribute wells

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            absorption (): Absorption of the species. Defaults to ListPlus()
            time (): Time of the measurement. Defaults to ListPlus()
            time_unit (): Unit of the time. Defaults to None
            reaction_volume (): Volume of the reaction. Defaults to None
            volume_unit (): Unit of the volume. Defaults to None
            init_conditions (): List of initial conditions of different species. Defaults to ListPlus()
            x_position (): X position of the well on the plate. Defaults to None
            y_position (): Y position of the well on the plate. Defaults to None
            wavelength (): Wavelength of the measurement. Defaults to None
        """

        params = {
            "absorption": absorption,
            "time": time,
            "time_unit": time_unit,
            "reaction_volume": reaction_volume,
            "volume_unit": volume_unit,
            "init_conditions": init_conditions,
            "x_position": x_position,
            "y_position": y_position,
            "wavelength": wavelength,
        }

        if id is not None:
            params["id"] = id

        self.wells.append(Well(**params))

        return self.wells[-1]
