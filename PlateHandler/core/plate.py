import sdRDM

from typing import List, Optional, Literal
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species
from .speciescondition import SpeciesCondition
from .speciestype import SpeciesType
from .well import Well
from PlateHandler.tools.spectramax_reader import read_spectramax


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

    species: List[Species] = Field(
        description="List of species present in wells of the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_wells(
        self,
        absorption: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        time_unit: Optional[str] = None,
        reaction_volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
        species_conditions: List[SpeciesCondition] = ListPlus(),
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
            species_conditions (): List of species conditions. Defaults to ListPlus()
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
            "species_conditions": species_conditions,
            "x_position": x_position,
            "y_position": y_position,
            "wavelength": wavelength,
        }

        if id is not None:
            params["id"] = id

        self.wells.append(Well(**params))

        return self.wells[-1]

    def add_to_species(
        self,
        type: SpeciesType,
        species_id: Optional[str] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            type (): Type of the species.
            species_id (): ID of the species. Defaults to None
            name (): Name of the species. Defaults to None
        """

        params = {
            "type": type,
            "species_id": species_id,
            "name": name,
        }

        if id is not None:
            params["id"] = id

        self.species.append(Species(**params))

        return self.species[-1]

    def assign_species(
            self,
            ids: List[str],
            species: Species,
            init_concs: List[float],
            conc_unit: str,
            to: Literal["all", "rows", "columns", "except"]
    ):

        cases = ["rows", "columns", "all", "except"]
        if not to in cases:
            raise AttributeError(
                f"Argument 'to' must be one of {cases}."
            )

        if to == "all":
            self.assign_species_to_all(species, init_concs[0], conc_unit)

        if to == "rows":
            self._species_to_rows(species, init_concs, conc_unit)

        return

    def assign_species_to_all(
        self,
        species: Species,
        init_conc: float,
        conc_unit: str
    ):
        nu_species = Species(**species.__dict__)
        nu_species.init_conc = init_conc
        nu_species.conc_unit = conc_unit

        for well in self.wells:
            well._update_species(nu_species)

    def assign_species_conditions_to_rows(
        self,
        row_ids: List[str],
        species: Species,
        init_concs: List[float],
        conc_unit: str,
    ):
        # Handle row_ids
        if not isinstance(row_ids, list):
            row_ids = [row_ids]

        if not all([isinstance(row_id, str) for row_id in row_ids]):
            raise AttributeError(
                "Argument 'row_ids' must be a list of strings."
            )

        # Handle init_concs
        if not isinstance(init_concs, list):
            init_concs = [init_concs]

        if len(init_concs) == 1:
            init_concs = init_concs * self.n_columns
        else:
            if len(init_concs) is not self.n_columns:
                raise AttributeError(
                    "Argument 'init_conc' must be a list of length equal to the number of columns."
                )

        for row_id in row_ids:
            for column_id, init_conc in zip(range(self.n_columns), init_concs):
                wells = self.get("wells", "id", f"{row_id}{column_id+1}")[0]
                for well in wells:
                    well.add_species_condition(
                        species_type=species,
                        init_conc=init_conc,
                        conc_unit=conc_unit,
                    )

    def assign_species_conditions_to_columns(self):
        raise NotImplementedError()

    @classmethod
    def from_file(
        cls,
        path: str,
        time: List[float],
        time_unit: str,
        ph: float = None,
        temperature: float = None,
        temperature_unit: str = None,
    ):
        return read_spectramax(cls, path, time, time_unit, ph, temperature, temperature_unit)
