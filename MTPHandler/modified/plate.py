import re
import sdRDM
import numpy as np

from typing import Dict, List, Optional, Literal, Union
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from MTPHandler.core import vessel
from MTPHandler.core.vessel import Vessel

from MTPHandler.modified.sboterm import SBOTerm


from .abstractspecies import AbstractSpecies
from .reactant import Reactant
from .protein import Protein
from .initcondition import InitCondition
from .well import Well
from .species import Species
from .speciestype import SpeciesType
from MTPHandler.tools.spectramax_reader import read_spectramax


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

    species: List[Union[AbstractSpecies, Protein, Reactant, None]] = Field(
        default_factory=ListPlus,
        multiple=True,
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

    def _add_to_species(self, new_species: AbstractSpecies) -> AbstractSpecies:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            type (): Type of the species.
            name (): Name of the species. Defaults to None
        """

        if any([species.id == new_species.id for species in self.species]):
            self.species = [new_species if species.id ==
                            new_species.id else species for species in self.species]

            return new_species

        else:
            self.species.append(new_species)

            return new_species

    def add_protein(
            self,
            id: str,
            name: str,
            constant: bool,
            sequence: str,
            **kwargs
    ):

        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "sequence": sequence,
            "vessel_id": vessel.id,
            **kwargs
        }

        return self._add_to_species(Protein(**params))

    def add_reactant(
            self,
            id: str,
            name: str,
            constant: bool,
            **kwargs
    ):

        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "vessel_id": vessel.id,
            **kwargs
        }

        return self._add_to_species(Reactant(**params))

    def _define_dummy_vessel(self):

        return Vessel(
            id="plate0",
            name="MTP 96 well",
            volume=200,
            unit="ul",
        )

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
        species: AbstractSpecies,
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

    def get_wells(
            self,
            ids: List[str] = None,
            rows: List[str] = None,
            columns: List[int] = None,
            wavelength: int = None,
    ) -> List[Well]:
        """
        Returns wells of a given wavelength, rows, columns, or individual ids
        of a plate.

        Args:
            ids (List[str], optional): Well ID (e.g. 'B10'). Defaults to None.
            rows (List[str], optional): ID of the rows (e.g. ['A', 'B']). Defaults to None.
            columns (List[int], optional): Column number. Starts at 1. Defaults to None.
            wavelength (int, optional): If plate was measured at different wavelengts,
            specify wavelength. Defaults to None.

        Raises:
            AttributeError: If multiple wavelengths were measured and no wavelength is specified.

        Returns:
            List[Well]: Wells matching to specified selection criteria.
        """

        # handel wavelength
        if not wavelength:
            if len(self.measured_wavelengths) == 1:
                wavelength = self.measured_wavelengths[0]
            else:
                raise AttributeError(
                    f"Argument 'wavelength' must be provided. Measured wavelengths are: {self.measured_wavelengths}"
                )

        # return wells, if ids are provided
        if ids and not any([rows, columns]):
            if not isinstance(ids, list):
                ids = [ids]

            return [self._get_well_by_id(id, wavelength) for id in ids]

        # return wells, if row_ids are provided
        elif rows and not any([ids, columns]):
            if not isinstance(rows, list):
                rows = [rows]

            well_rows = [self._get_wells_by_row_id(
                row_id, wavelength) for row_id in rows]

            return [well for row in well_rows for well in row]

        # return wells, if column_ids are provided
        elif columns and not any([ids, rows]):
            if not isinstance(columns, list):
                columns = [columns]

            well_columns = [self._get_wells_by_column_id(
                column_id, wavelength) for column_id in columns]

            return [well for column in well_columns for well in column]

        # return all wells of a given wavelength
        else:
            return [well for well in self.wells if well.wavelength == wavelength]

    def _get_well_by_id(self, id: str, wavelength: int) -> Well:

        for well in self.wells:
            if well.id == id and well.wavelength == wavelength:
                return well

        raise ValueError(f"No well found with id {id}")

    def _get_wells_by_column_id(self, column_id: int, wavelength: int) -> Well:

        x_position = column_id-1
        y_positions = [well.y_position for well in self.wells if well.x_position ==
                       x_position and well.wavelength == wavelength]

        return [self._get_well_by_xy(x_position, y_pos, wavelength) for y_pos in y_positions]

    def _get_wells_by_row_id(self, row_id: str, wavelength: int) -> List[Well]:

        y_position = ord(row_id)-65
        x_positions = [well.x_position for well in self.wells if well.y_position ==
                       y_position and well.wavelength == wavelength]

        return [self._get_well_by_xy(x_pos, y_position, wavelength) for x_pos in x_positions]

    def _get_well_by_xy(self, x_position: int, y_position: int, wavelength: int) -> Well:

        for well in self.wells:
            if well.x_position == x_position and well.y_position == y_position and well.wavelength == wavelength:
                return well

        raise ValueError(
            f"No well found with x position {x_position} and y position {y_position}")

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

        # validate well_id
        self._validate_well_id(well_ids)

        wells = (well for well in self.wells if well.id not in well_ids)
        for well in wells:
            well.add_species_condition(
                species_type=species,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

    def get_well(self, _id: str) -> Well:

        for well in self.wells:
            if well.id == _id:
                return well

        raise ValueError(f"No well found with id {_id}")

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
        wells = self.get_wells(wavelength=wavelength)

        # get wells for blanking
        blank_wells = self._get_blanks(wells, species)

        # calculate mean absorption across all measured blank values
        species_abso_contribution = [
            np.mean([well.absorption for well in blank_wells])][0]

        # apply to wells where species is present
        for well in wells:
            if any(condition.species_id == species.id for condition in well.init_conditions):
                well.absorption = [
                    absorption - species_abso_contribution for absorption in well.absorption]
                blanked_species = well._get_species_condition(species)
                blanked_species.was_blanked = True

    def get_species(self, _id: str) -> Species:

        for species in self.species:
            if species.id == _id:
                return species

        raise ValueError(f"No species found with id {_id}")

    def _get_catalyst(self) -> bool:

        for species in self.species:
            if species.ontology == SBOTerm.CATALYST.value:
                return species

        return None

    @staticmethod
    def _get_blanks(wells: List[Well], species: Species) -> List[Well]:

        blank_wells = []

        for well in wells:
            species_blanked_status = [
                condition.was_blanked for condition in well.init_conditions]
            if not species_blanked_status.count(False) == 1:
                continue

            if well.init_conditions[species_blanked_status.index(False)].species_id == species.id:
                blank_wells.append(well)

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
