import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from CaliPytion import Standard

from .abstractspecies import AbstractSpecies
from .initcondition import InitCondition
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
        species: Optional[AbstractSpecies] = None,
        init_conc: Optional[float] = None,
        conc_unit: Optional[str] = None,
        was_blanked: bool = False,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'InitCondition' to attribute init_conditions

        Args:
            id (str): Unique identifier of the 'InitCondition' object. Defaults to 'None'.
            species (): Reference to species. Defaults to None
            init_conc (): Initial concentration of the species. Defaults to None
            conc_unit (): Concentration unit. Defaults to None
            was_blanked (): Whether the species' absorption contribution was subtracted from the absorption signal. Defaults to False
        """

        params = {
            "species_id": species.id,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
            "was_blanked": was_blanked,
        }

        new_condition = InitCondition(**params)

        print(f"assigning to {self.id} {init_conc}")

        if any([condition.species_id == new_condition.species_id for condition in self.init_conditions]):
            self.init_conditions = [new_condition if condition.species_id ==
                                    new_condition.species_id else condition for condition in self.init_conditions]

        else:
            self.init_conditions.append(new_condition)

    def _contains_species(self, species: AbstractSpecies) -> bool:

        for condition in self.init_conditions:
            if condition.species_id == species.id:
                return True

        return False

    def _is_blanked(self, species_id: str) -> bool:

        other_species_blanked = []
        for condition in self.init_conditions:
            if condition.species_id == species_id:
                species_was_blanked = condition.was_blanked

            else:
                other_species_blanked.append(condition.was_blanked)

        if not species_was_blanked and all(other_species_blanked):
            return True
        else:
            return False

    def _get_species_condition(self, species: AbstractSpecies) -> InitCondition:

        for condition in self.init_conditions:
            if condition.species_id == species.id:
                return condition

        raise ValueError(f"Species {species} not found in well {self.id}")

    def to_concentration(self, standard: Standard, **kwargs) -> list[float]:

        if not standard.model_result.was_fitted:
            raise ValueError(
                f"Standard {standard.id} was not fitted for species {standard.species_id}")

        if not self._is_blanked(standard.species_id):
            raise ValueError(
                f"Well {self.id} was not blanked for species {standard.species_id}")

        if not self.wavelength == standard.wavelength:
            raise ValueError(
                f"Standard at {standard.wavelength} nm not applicable for well {self.id} measured at {self.wavelength}")

        return standard.model_result.calculate(self.absorption, **kwargs)
