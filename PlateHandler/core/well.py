import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species
from .speciescondition import SpeciesCondition
from .speciestype import SpeciesType


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

    species_conditions: List[SpeciesCondition] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="List of species conditions",
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

    def __init__(self, **data):
        super().__init__(**data)
        self._add_buffer()

    def add_species_condition(
        self,
        species_type: Optional[Species] = None,
        init_conc: Optional[float] = None,
        conc_unit: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'SpeciesCondition' to attribute species_conditions

        Args:
            id (str): Unique identifier of the 'SpeciesCondition' object. Defaults to 'None'.
            species_type (): Reference to species. Defaults to None
            init_conc (): Initial concentration of the species. Defaults to None
            conc_unit (): Concentration unit. Defaults to None
        """

        params = {
            "species_type": species_type,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
        }

        new_condition = SpeciesCondition(**params)

        if any([condition.species_type == new_condition.species_type for condition in self.species_conditions]):
            self.species_conditions = [new_condition if condition.species_type ==
                                       new_condition.species_type else condition for condition in self.species_conditions]

        else:
            self.species_conditions.append(new_condition)

    def _add_buffer(self):
        '''Adds a species with `SpeciesType.BUFFER` to the well'''

        self.add_species_condition(
            species_type=SpeciesType.BUFFER
        )
