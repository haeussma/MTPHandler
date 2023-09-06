import sdRDM

from typing import Optional, Union, List
from pydantic import Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species
from .well import Well


@forge_signature
class Blank(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("blankINDEX"),
        xml="@id",
    )

    well_ids: List[Union[Well, str]] = Field(
        reference="Well.id",
        default_factory=ListPlus,
        multiple=True,
        description="ID of the well",
    )

    species_id: Union[Species, str, None] = Field(
        default=Species(),
        reference="Species.id",
        description="ID of the species",
    )

    concentrations: List[float] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Concentration of the species",
    )

    concentration_unit: Optional[str] = Field(
        default=None,
        description="Unit of the concentration",
    )

    absorption_contribution: Optional[float] = Field(
        default=None,
        description="Absorption contribution of the species",
    )

    def add_to_well_ids(
        self,
        absorption: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        time_unit: Optional[str] = None,
        reaction_volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
        species: Optional[Species] = None,
        x_position: Optional[int] = None,
        y_position: Optional[int] = None,
        species_id: Optional[str] = None,
        wavelength: Optional[int] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Well' to attribute well_ids

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            absorption (): Absorption of the species. Defaults to ListPlus()
            time (): Time of the measurement. Defaults to ListPlus()
            time_unit (): Unit of the time. Defaults to None
            reaction_volume (): Volume of the reaction. Defaults to None
            volume_unit (): Unit of the volume. Defaults to None
            species (): Species present in the well.. Defaults to None
            x_position (): X position of the well on the plate. Defaults to None
            y_position (): Y position of the well on the plate. Defaults to None
            species_id (): ID of the species. Defaults to None
            wavelength (): Wavelength of the measurement. Defaults to None
        """

        params = {
            "absorption": absorption,
            "time": time,
            "time_unit": time_unit,
            "reaction_volume": reaction_volume,
            "volume_unit": volume_unit,
            "species": species,
            "x_position": x_position,
            "y_position": y_position,
            "species_id": species_id,
            "wavelength": wavelength,
        }

        if id is not None:
            params["id"] = id

        self.well_ids.append(Well(**params))

        return self.well_ids[-1]

    @validator("well_ids")
    def get_well_ids_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""

        from .well import Well

        if isinstance(value, Well):
            return value.id
        elif isinstance(value, str):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [Well, str] got '{type(value).__name__}' instead."
            )

    @validator("species_id")
    def get_species_id_reference(cls, value):
        """Extracts the ID from a given object to create a reference"""

        from .species import Species

        if isinstance(value, Species):
            return value.id
        elif isinstance(value, str):
            return value
        elif value is None:
            return value
        else:
            raise TypeError(
                f"Expected types [Species, str] got '{type(value).__name__}' instead."
            )
