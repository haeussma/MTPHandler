import sdRDM

import itertools as it
import re
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional, Union, List, Callable, Dict, Literal
from pydantic import Field, StrictBool
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from datetime import datetime as Datetime
from collections import defaultdict
from types import NoneType
from plotly.subplots import make_subplots
from CaliPytion import Calibrator, Standard
from MTPHandler.ioutils import initialize_calibrator, create_enzymeml
from .abstractspecies import AbstractSpecies
from .initcondition import InitCondition
from .well import Well
from .photometricmeasurement import PhotometricMeasurement
from .vessel import Vessel
from .protein import Protein
from .reactant import Reactant
from .sboterm import SBOTerm


@forge_signature
class Plate(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("plateINDEX"),
        xml="@id",
    )

    n_rows: int = Field(..., description="Number of rows on the plate")

    n_columns: int = Field(..., description="Number of columns on the plate")

    date_measured: Optional[Datetime] = Field(
        default=None, description="Date and time when the plate was measured"
    )

    times: List[float] = Field(
        description=(
            "Time points of the measurement, corresponding to temperature measurements"
        ),
        default_factory=ListPlus,
        multiple=True,
    )

    time_unit: Optional[str] = Field(default=None, description="Unit of the time")

    temperatures: List[float] = Field(
        description="Thermostat temperature", multiple=True, default_factory=ListPlus
    )

    temperature_unit: str = Field(..., description="Unit of the temperature")

    max_volume: Optional[float] = Field(
        default=None, description="Maximum volume of the wells"
    )

    max_volume_unit: Optional[str] = Field(
        default=None, description="Unit of the maximum volume"
    )

    ph: float = Field(..., description="pH of the reaction")

    wells: List[Well] = Field(
        description="List of wells on the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    measured_wavelengths: List[float] = Field(
        description="Measured wavelengths", default_factory=ListPlus, multiple=True
    )

    wavelength_unit: Optional[str] = Field(
        default=None, description="Unit of the wavelength"
    )

    species: List[AbstractSpecies] = Field(
        description="List of species present in wells of the plate",
        default_factory=ListPlus,
        multiple=True,
    )

    def add_to_wells(
        self,
        x_position: int,
        y_position: int,
        init_conditions: List[InitCondition] = ListPlus(),
        measurements: List[PhotometricMeasurement] = ListPlus(),
        volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Well' to attribute wells

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            x_position (): X position of the well on the plate.
            y_position (): Y position of the well on the plate.
            init_conditions (): List of initial conditions of different species. Defaults to ListPlus()
            measurements (): List of photometric measurements. Defaults to ListPlus()
            volume (): Volume of the reaction. Defaults to None
            volume_unit (): Unit of the volume. Defaults to None
        """
        params = {
            "x_position": x_position,
            "y_position": y_position,
            "init_conditions": init_conditions,
            "measurements": measurements,
            "volume": volume,
            "volume_unit": volume_unit,
        }
        if id is not None:
            params["id"] = id
        self.wells.append(Well(**params))
        return self.wells[-1]

    def add_to_species(
        self,
        name: str,
        vessel_id: Vessel,
        constant: StrictBool,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'AbstractSpecies' to attribute species

        Args:
            id (str): Unique identifier of the 'AbstractSpecies' object. Defaults to 'None'.
            name (): None.
            vessel_id (): None.
            constant (): None.
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """
        params = {
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }
        if id is not None:
            params["id"] = id
        self.species.append(AbstractSpecies(**params))
        return self.species[-1]

    def add_abstract_species_to_species(
        self,
        name: str,
        vessel_id: Vessel,
        constant: StrictBool,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'AbstractSpecies' to attribute species

        Args:
            id (str): Unique identifier of the 'AbstractSpecies' object. Defaults to 'None'.
            name (): None.
            vessel_id (): None.
            constant (): None.
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """
        params = {
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }
        if id is not None:
            params["id"] = id
        self.species.append(AbstractSpecies(**params))
        return self.species[-1]

    def add_protein_to_species(
        self,
        sequence: str,
        name: str,
        vessel_id: Vessel,
        constant: StrictBool,
        ecnumber: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[str] = None,
        uniprotid: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.CATALYST,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Protein' to attribute species

        Args:
            id (str): Unique identifier of the 'Protein' object. Defaults to 'None'.
            sequence (): Amino acid sequence of the protein.
            name (): None.
            vessel_id (): None.
            constant (): None.
            ecnumber (): EC number of the protein.. Defaults to None
            organism (): Organism the protein was expressed in.. Defaults to None
            organism_tax_id (): Taxonomy identifier of the expression host.. Defaults to None
            uniprotid (): Unique identifier referencing a protein entry at UniProt. Use this identifier to initialize the object from the UniProt database.. Defaults to None
            ontology (): None. Defaults to SBOTerm.CATALYST
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """
        params = {
            "sequence": sequence,
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "ecnumber": ecnumber,
            "organism": organism,
            "organism_tax_id": organism_tax_id,
            "uniprotid": uniprotid,
            "ontology": ontology,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }
        if id is not None:
            params["id"] = id
        self.species.append(Protein(**params))
        return self.species[-1]

    def add_reactant_to_species(
        self,
        name: str,
        vessel_id: Vessel,
        constant: StrictBool,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        chebi_id: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.SMALL_MOLECULE,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Reactant' to attribute species

        Args:
            id (str): Unique identifier of the 'Reactant' object. Defaults to 'None'.
            name (): None.
            vessel_id (): None.
            constant (): None.
            smiles (): Simplified Molecular Input Line Entry System (SMILES) encoding of the reactant.. Defaults to None
            inchi (): International Chemical Identifier (InChI) encoding of the reactant.. Defaults to None
            chebi_id (): Unique identifier of the CHEBI database. Use this identifier to initialize the object from the CHEBI database.. Defaults to None
            ontology (): None. Defaults to SBOTerm.SMALL_MOLECULE
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """
        params = {
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "smiles": smiles,
            "inchi": inchi,
            "chebi_id": chebi_id,
            "ontology": ontology,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }
        if id is not None:
            params["id"] = id
        self.species.append(Reactant(**params))
        return self.species[-1]

    def add_protein(self, id: str, name: str, constant: bool, sequence: str, **kwargs):
        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "sequence": sequence,
            "vessel_id": vessel.id,
            **kwargs,
        }

        return self._add_species(Protein(**params))

    def _add_species(self, new_species: AbstractSpecies) -> AbstractSpecies:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            type (): Type of the species.
            name (): Name of the species. Defaults to None
        """

        if any([species.id == new_species.id for species in self.species]):
            self.species = [
                new_species if species.id == new_species.id else species
                for species in self.species
            ]

            return new_species

        else:
            self.species.append(new_species)

            return new_species

    def add_reactant(self, id: str, name: str, constant: bool, **kwargs):
        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "vessel_id": vessel.id,
            **kwargs,
        }

        print(params)

        return self._add_species(Reactant(**params))

    def _define_dummy_vessel(self):
        return Vessel(
            id="plate0",
            name="MTP 96 well",
            volume=200,
            unit="ul",
        )

    def assign_species(
        self,
        species: AbstractSpecies,
        init_conc: Union[float, List[float]],
        conc_unit: str,
        to: Literal["all", "rows", "columns", "except"],
        ids: Union[str, List[str], int, List[int]] = None,
    ):
        cases = ["rows", "columns", "all", "except"]
        if not to in cases:
            raise AttributeError(f"Argument 'to' must be one of {cases}.")

        if not isinstance(init_conc, list):
            init_conc = [init_conc]

        if not isinstance(ids, list) and not isinstance(ids, NoneType):
            ids = [ids]

        if to == "all":
            self.assign_species_to_all(species, init_conc, conc_unit)

        elif to == "columns":
            self.assign_species_to_columns(
                column_ids=ids,
                species=species,
                init_concs=init_conc,
                conc_unit=conc_unit,
            )

        elif to == "rows":
            self.assign_species_to_rows(
                row_ids=ids, species=species, init_concs=init_conc, conc_unit=conc_unit
            )

        else:
            self.assign_species_to_all_except(
                well_ids=ids, species=species, init_conc=init_conc, conc_unit=conc_unit
            )

        return

    def assign_species_to_all(
        self,
        species: AbstractSpecies,
        init_conc: float,
        conc_unit: str,
    ):
        if not len(init_conc) == 1:
            raise AttributeError(
                "Argument 'init_conc' must be a float, when assigning to all wells."
            )

        for well in self.wells:
            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc[0],
                conc_unit=conc_unit,
            )
        print(f"Assigned {species.name} to all wells.")

    def assign_species_to_columns(
        self,
        column_ids: List[int],
        species: AbstractSpecies,
        init_concs: List[float],
        conc_unit: str,
    ):
        # Handle column_ids
        if not all([isinstance(column_id, int) for column_id in column_ids]):
            raise AttributeError("Argument 'column_ids' must be a list of integers.")

        if not all([column_id <= self.n_columns for column_id in column_ids]):
            raise AttributeError(
                "Argument 'column_ids' must be a list of integers between 1 and"
                f" {self.n_columns+1}."
            )

        # Handle init_concs
        if len(init_concs) == 1:
            init_concs = init_concs * self.n_rows

        if not len(init_concs) == self.n_rows:
            raise AttributeError(
                f"Argument 'init_concs' must be a list of length {self.n_rows}."
            )

        for column_id in column_ids:
            for row_id, init_conc in zip(range(self.n_rows), init_concs):
                [
                    well.add_to_init_conditions(
                        species_id=species.id,
                        init_conc=init_conc,
                        conc_unit=conc_unit,
                    )
                    for well in self.wells
                    if well.y_position == row_id and well.x_position == column_id - 1
                ]

        print(
            f"Assigned {species.name} with concentrations of"
            f" {init_concs} {conc_unit} to columns {column_ids}."
        )

    def assign_species_to_rows(
        self,
        row_ids: List[str],
        species: AbstractSpecies,
        init_concs: List[float],
        conc_unit: str,
    ):
        # Handle row_ids
        if not row_ids:
            row_ids = [chr(i) for i in range(65, 65 + self.n_rows)]

        if not all([isinstance(row_id, str) for row_id in row_ids]):
            raise AttributeError("Argument 'row_ids' must be a list of strings.")

        else:
            row_ids = [_id.upper() for _id in row_ids]

        # Handle init_concs
        if len(init_concs) == 1:
            init_concs = init_concs * self.n_columns

        if not len(init_concs) == self.n_columns:
            raise AttributeError(
                f"Argument 'init_concs' must be a list of length {self.n_columns}."
            )

        # Add concentration array to wells in rows
        for row_id in row_ids:
            for column_id, init_conc in zip(range(self.n_columns), init_concs):
                [
                    well.add_to_init_conditions(
                        species_id=species.id,
                        init_conc=init_conc,
                        conc_unit=conc_unit,
                    )
                    for well in self.wells
                    if well.x_position == column_id
                    and well.y_position == ord(row_id) - 65
                ]

        print(
            f"Assigned {species.name} with concentrations of"
            f" {init_concs} {conc_unit} to rows {row_ids}."
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
                    "Argument 'wavelength' must be provided. Measured wavelengths are:"
                    f" {self.measured_wavelengths}"
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

            well_rows = [
                self._get_wells_by_row_id(row_id, wavelength) for row_id in rows
            ]

            return [well for row in well_rows for well in row]

        # return wells, if column_ids are provided
        elif columns and not any([ids, rows]):
            if not isinstance(columns, list):
                columns = [columns]

            well_columns = [
                self._get_wells_by_column_id(column_id, wavelength)
                for column_id in columns
            ]

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
        x_position = column_id - 1
        y_positions = [
            well.y_position
            for well in self.wells
            if well.x_position == x_position and well.wavelength == wavelength
        ]

        return [
            self._get_well_by_xy(x_position, y_pos, wavelength) for y_pos in y_positions
        ]

    def _get_wells_by_row_id(self, row_id: str, wavelength: int) -> List[Well]:
        y_position = ord(row_id) - 65
        x_positions = [
            well.x_position
            for well in self.wells
            if well.y_position == y_position and well.wavelength == wavelength
        ]

        return [
            self._get_well_by_xy(x_pos, y_position, wavelength) for x_pos in x_positions
        ]

    def _get_well_by_xy(
        self, x_position: int, y_position: int, wavelength: int
    ) -> Well:
        for well in self.wells:
            if (
                well.x_position == x_position
                and well.y_position == y_position
                and well.wavelength == wavelength
            ):
                return well

        raise ValueError(
            f"No well found with x position {x_position} and y position {y_position}"
        )

    def assign_species_to_all_except(
        self,
        well_ids: List[str],
        species: AbstractSpecies,
        init_conc: float,
        conc_unit: str,
    ):
        if not len(init_conc) == 1:
            raise AttributeError(
                "Argument 'init_conc' must be a float, when assigning to all wells."
            )

        # validate well_id
        self._validate_well_id(well_ids)

        wells = (well for well in self.wells if well.id not in well_ids)
        for well in wells:
            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc[0],
                conc_unit=conc_unit,
            )

    def get_well(self, _id: str, wavelength: float) -> Well:
        for well in self.wells:
            if well.id == _id and well.wavelength == wavelength:
                return well

        raise ValueError(f"No well found with id {_id} and wavelength {wavelength}")

    def calibrate(
        self,
        species: AbstractSpecies,
        wavelength: int,
        cutoff: float = None,
    ) -> Calibrator:
        return initialize_calibrator(
            plate=self,
            species=species,
            wavelength=wavelength,
            cutoff=cutoff,
        )

    def to_enzymeml(
        self,
        name: str,
        detected_reactant: Reactant,
        mapping_reactant: Reactant = None,
        reactant_standard: Standard = None,
        wavelength: int = None,
        path: str = None,
    ) -> "EnzymeML":
        return create_enzymeml(
            name=name,
            plate=self,
            detected_reactant=detected_reactant,
            mapping_reactant=mapping_reactant,
            reactant_standard=reactant_standard,
            wavelength=wavelength,
            path=path,
        )

    def _already_blanked(self, wells: List[Well], species: AbstractSpecies) -> bool:
        wells = [well for well in wells if well._contains_species(species.id)]

        wells = [
            well
            for well in wells
            if well._get_species_condition(species.id).init_conc != 0
        ]

        if all([well._get_species_condition(species.id).was_blanked for well in wells]):
            return True
        else:
            return False

    def blank_species(
        self,
        species: AbstractSpecies,
        wavelength: int,
    ):
        # get wells with specified wavelength:
        if not wavelength in self.measured_wavelengths:
            raise AttributeError(
                f"Plate does not contain measurements with wavelength {wavelength} nm."
            )
        wells = self.get_wells(wavelength=wavelength)

        if self._already_blanked(wells=wells, species=species):
            print(f"{species.name} was already blanked.")
            return

        # get wells for blanking
        blank_wells = self._get_blanks(wells=wells, species=species)

        # get mapping of concentration to blank wells
        blank_well_mapping = self._get_blank_conc_mapping(
            wells=blank_wells, species=species
        )

        conc_mean_blank_mapping = {}
        for conc, blanking_wells in blank_well_mapping.items():
            mean_absorption = np.nanmean([well.absorption for well in blanking_wells])
            conc_mean_blank_mapping[conc] = mean_absorption

        # apply to wells where species is present in respective concentration
        blanked_wells = []
        for well in wells:
            if species.id not in [
                condition.species_id for condition in well.init_conditions
            ]:
                continue

            condition = well._get_species_condition(species.id)

            if condition.was_blanked:
                continue

            if condition.init_conc not in conc_mean_blank_mapping.keys():
                print(
                    f"Well {well.id} was not blanked for initial"
                    f" {species.name} concentration"
                    f" {condition.init_conc} ({condition.conc_unit})."
                )
                continue

            prior = well.absorption[3]

            well.absorption = [
                absorption - conc_mean_blank_mapping[condition.init_conc]
                for absorption in well.absorption
            ]

            blanked_wells.append(well.id)

            condition.was_blanked = True

        print(f"Blanked {len(blanked_wells)} wells containing {species.name}.")

    def visualize(self, zoom: bool = False):
        if zoom:
            shared_yaxes = False
        else:
            shared_yaxes = True

        fig = make_subplots(
            rows=self.n_rows,
            cols=self.n_columns,
            shared_xaxes=True,
            subplot_titles=[well.id for well in self.wells],
            shared_yaxes=shared_yaxes,
        )
        colors = px.colors.qualitative.Plotly

        for well in self.wells:
            for measurement, color in zip(well.measurements, colors):
                fig.add_trace(
                    go.Scatter(
                        x=self.times,
                        y=measurement.absorptions,
                        name=f"{measurement.wavelength} nm",
                        showlegend=False,
                        line=dict(color=color),
                    ),
                    col=well.x_position + 1,
                    row=well.y_position + 1,
                )

        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)
        fig.update_traces(hovertemplate="%{y:.2f}")

        fig.update_layout(
            plot_bgcolor="white",
            hovermode="x",
            title=dict(
                text=f"pH {self.ph}, {self.temperatures[0]} °{self.temperature_unit}",
            ),
            margin=dict(l=20, r=20, t=100, b=20),
        )

        fig.show()

    def get_species(self, _id: str) -> AbstractSpecies:
        for species in self.species:
            if species.id == _id:
                return species

        raise ValueError(f"No species found with id {_id}")

    def _generate_possible_well_ids(self):
        characters = "ABCDEFGH"
        integers = range(1, 13)  # 1 to 12

        sub_char = characters[: self.n_rows]
        sub_int = integers[: self.n_columns]

        # Generate combinations of characters and integers
        combinations = [
            "".join(item) for item in it.product(sub_char, map(str, sub_int))
        ]

        return combinations

    @staticmethod
    def _get_blanks(wells: List[Well], species: AbstractSpecies) -> List[Well]:
        blank_wells = []
        for well in wells:
            # ignore virtual species with concentration of 0
            conditions = [
                condition
                for condition in well.init_conditions
                if condition.init_conc != 0
            ]

            # ignore species that were already blanked
            conditions = [
                condition for condition in conditions if condition.was_blanked == False
            ]

            # if not exactly one species was not blanked, continue
            if not len(conditions) == 1:
                continue

            # if the species is not the one to be blanked, continue
            if not conditions[0].species_id == species.id:
                continue

            blank_wells.append(well)

        if len(blank_wells) == 0:
            raise AttributeError(
                "No wells for calculating the blank found for"
                f" {species.name} ({species.id}). You might need to blank another"
                " species first."
            )

        return blank_wells

    def _get_blank_conc_mapping(
        self, wells: List[Well], species: AbstractSpecies
    ) -> Dict[float, List[Well]]:
        blank_conc_mapping = defaultdict(list)
        for well in wells:
            condition = well._get_species_condition(species.id)

            blank_conc_mapping[condition.init_conc].append(well)

        return blank_conc_mapping

    @staticmethod
    def _validate_well_id(well_ids) -> bool:
        WELL_ID = re.compile(r"[A-Z][1-9]\d?")

        if not all(WELL_ID.match(well_id) for well_id in well_ids):
            raise ValueError(
                "Invalid well id(s) provided:"
                f" {[well for well in well_ids if not WELL_ID.match(well)]}"
            )

    @classmethod
    def from_reader(cls, reader: Callable, path: str, **kwargs):
        return reader(cls, path, **kwargs)
