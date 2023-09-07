import re
import sdRDM
import numpy as np

from typing import Dict, List, Optional, Literal
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .initcondition import InitCondition
from .well import Well
from .species import Species
from .speciestype import SpeciesType
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

    measured_wavelengths: List[int] = Field(
        description="Measured wavelengths in nm",
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

    def add_to_species(
        self, type: SpeciesType, name: Optional[str] = None, id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            type (): Type of the species.
            name (): Name of the species. Defaults to None
        """

        params = {
            "type": type,
            "name": name,
            "id": id,
        }

        new_species = Species(**params)

        if any([species.id == new_species.id for species in self.species]):
            self.species = [new_species if species.id ==
                            new_species.id else species for species in self.species]

            return new_species

        else:
            self.species.append(new_species)

            return new_species

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

                # if the concentration of a species is 0, it is not added, since its not present
                if init_conc == 0:
                    continue

                wells = (well for well in self.wells if well.id ==
                         f"{row_id}{column_id+1}")
                for well in wells:
                    well.add_to_init_conditions(
                        species=species,
                        init_conc=init_conc,
                        conc_unit=conc_unit,
                    )

    def assign_species_conditions_to_all_except(
            self,
            well_ids: List[str],
            species: Species,
            init_conc: float,
            conc_unit: str
    ):

        if not isinstance(well_ids, list):
            well_ids = [well_ids]

        if not isinstance(init_conc, float):
            raise AttributeError(
                "Argument 'init_conc' must be a float."
            )

        # Check well_id correctness
        self._validate_well_id(well_ids)

        wells = (well for well in self.wells if well.id not in well_ids)
        for well in wells:
            well.add_species_condition(
                species_type=species,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

    def blank_species(
            self,
            species: Species,
            wavelength: int,
    ):

        # get wells with specified wavelength:
        if not wavelength in self.measured_wavelengths:
            raise AttributeError(
                f"Plate does not contain measurements with wavelength {wavelength} nm."
            )
        wells = [well for well in self.wells if well.wavelength == wavelength]

        # get wells for blanking
        blank_wells = self._get_blanks_for_species(wells, species)

        # calculate mean absorption across all measured blank values
        species_abso_contribution = [
            np.mean([well.absorption for well in blank_wells])][0]

        # apply to wells where species is present
        for well in wells:
            if any(condition.species == species.id for condition in well.init_conditions):
                well.absorption = [
                    absorption - species_abso_contribution for absorption in well.absorption]
                blanked_species = well._get_species_condition(species)
                blanked_species.was_blanked = True

    @staticmethod
    def _get_blanks_for_species(wells: List[Well], species: Species) -> List[Well]:

        blank_wells = []

        for well in wells:
            species_blanked_status = [
                condition.was_blanked for condition in well.init_conditions]
            if not species_blanked_status.count(False) == 1:
                continue

            if well.init_conditions[species_blanked_status.index(False)].species == species.id:
                blank_wells.append(well)
                print(
                    f"added {well.id}, with condition length: {len(well.init_conditions)} and [0] being {well.init_conditions[species_blanked_status.index(False)].species}")

        return blank_wells

    @staticmethod
    def _validate_well_id(well_ids) -> bool:

        WELL_ID = re.compile(r"[A-Z][1-9]\d?")

        if not all(WELL_ID.match(well_id) for well_id in well_ids):
            raise ValueError(
                f"Invalid well id(s) provided: {[well for well in well_ids if not WELL_ID.match(well)]}"
            )

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
