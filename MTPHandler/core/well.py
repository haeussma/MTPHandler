import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .initcondition import InitCondition
from .abstractspecies import AbstractSpecies


@forge_signature
class Well(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("wellINDEX"),
        xml="@id",
    )

    absorption: List[float] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Absorption of the species",
    )

    time: List[float] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Time of the measurement",
    )

    time_unit: Optional[str] = Field(
        default=None,
        description="Unit of the time",
    )

    reaction_volume: Optional[float] = Field(
        default=None,
        description="Volume of the reaction",
    )

    volume_unit: Optional[str] = Field(
        default=None,
        description="Unit of the volume",
    )

    init_conditions: List[InitCondition] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="List of initial conditions of different species",
    )

    x_position: Optional[int] = Field(
        default=None,
        description="X position of the well on the plate",
    )

    y_position: Optional[int] = Field(
        default=None,
        description="Y position of the well on the plate",
    )

    wavelength: Optional[int] = Field(
        default=None,
        description="Wavelength of the measurement",
    )

    def add_to_init_conditions(
        self,
        species_id: Optional[AbstractSpecies] = None,
        init_conc: Optional[float] = None,
        conc_unit: Optional[str] = None,
        was_blanked: bool = False,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'InitCondition' to attribute init_conditions

        Args:
            id (str): Unique identifier of the 'InitCondition' object. Defaults to 'None'.
            species_id (): Reference to species. Defaults to None
            init_conc (): Initial concentration of the species. Defaults to None
            conc_unit (): Concentration unit. Defaults to None
            was_blanked (): Whether the species' absorption contribution was subtracted from the absorption signal. Defaults to False
        """

        params = {
            "species_id": species_id,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
            "was_blanked": was_blanked,
        }

        if id is not None:
            params["id"] = id

        self.init_conditions.append(InitCondition(**params))

        return self.init_conditions[-1]
