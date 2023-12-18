import sdRDM

import itertools as it
import re
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional, Union, List, Callable, Dict, Literal
from pydantic import Field, StrictBool
from pydantic import StrictBool, Field
from pydantic import Field, StrictBool
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from datetime import datetime as Datetime
from collections import defaultdict
from types import NoneType
from plotly.subplots import make_subplots
from CaliPytion.core import Calibrator, Standard
from MTPHandler.ioutils import initialize_calibrator
from MTPHandler.ioutils.enzymeml import sort_measurements
from .well import Well
from .photometricmeasurement import PhotometricMeasurement
from .initcondition import InitCondition
from .abstractspecies import AbstractSpecies
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

        return self._add_species(Reactant(**params))

    def _define_dummy_vessel(self):
        return Vessel(
            id="v0",
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
        contributes_to_signal: Optional[bool] = None,
    ):
        cases = ["rows", "columns", "all", "except"]
        if not to in cases:
            raise AttributeError(f"Argument 'to' must be one of {cases}.")

        if not isinstance(init_conc, list):
            init_conc = [init_conc]

        if not isinstance(ids, list) and not isinstance(ids, NoneType):
            ids = [ids]

        if to == "all":
            self.assign_species_to_all(
                species=species,
                init_conc=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        elif to == "columns":
            self.assign_species_to_columns(
                column_ids=ids,
                species=species,
                init_concs=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        elif to == "rows":
            self.assign_species_to_rows(
                row_ids=ids,
                species=species,
                init_concs=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        else:
            self.assign_species_to_all_except(
                well_ids=ids,
                species=species,
                init_conc=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        return

    def assign_species_to_all(
        self,
        species: AbstractSpecies,
        init_conc: float,
        conc_unit: str,
        contributes_to_signal: Optional[bool] = None,
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

            if contributes_to_signal is not None:
                contributes = contributes_to_signal
            else:
                if init_conc == 0:
                    contributes = False
                else:
                    contributes = True

            for measurement in well.measurements:
                measurement.add_to_blank_states(
                    species_id=species.id,
                    contributes_to_signal=contributes,
                )

        print(f"Assigned {species.name} to all wells.")

    def assign_species_to_columns(
        self,
        column_ids: List[int],
        species: AbstractSpecies,
        init_concs: List[float],
        conc_unit: str,
        contributes_to_signal: Optional[bool] = None,
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
                for well in self.wells:
                    if not well.x_position == column_id - 1:
                        continue
                    if not well.y_position == row_id:
                        continue

                    well.add_to_init_conditions(
                        species_id=species.id, init_conc=init_conc, conc_unit=conc_unit
                    )

                    if contributes_to_signal is not None:
                        contributes = contributes_to_signal
                    else:
                        if init_conc == 0:
                            contributes = False
                        else:
                            contributes = True

                    for measurement in well.measurements:
                        measurement.add_to_blank_states(
                            species_id=species.id,
                            contributes_to_signal=contributes,
                        )

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
        contributes_to_signal: Optional[bool] = None,
    ):
        # Handle row_ids
        if not row_ids:
            row_ids = [chr(i) for i in range(65, 65 + self.n_rows)]

        if not all([isinstance(row_id, str) for row_id in row_ids]):
            raise AttributeError("Argument 'row_ids' must be a list of strings.")

        else:
            row_ids = [_id.upper() for _id in row_ids]

        # Handle init_concs

        for row_id in row_ids:
            if len(init_concs) == 1:
                init_concs = init_concs * len(self._get_wells_by_row_id(row_id))

            self._assign_species_to_row(
                row_id=row_id,
                init_concs=init_concs,
                species=species,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

    def _assign_species_to_row(
        self,
        row_id: str,
        init_concs: List[float],
        species: AbstractSpecies,
        conc_unit: str,
        contributes_to_signal: Optional[bool],
    ):
        wells = self._get_wells_by_row_id(row_id)

        if not len(init_concs) == len(wells):
            raise AttributeError(
                f"Number of initial concentrations ({len(init_concs)}) does not match"
                f"number of rows in row {row_id} ({len(wells)})"
            )

        for well, init_conc in zip(wells, init_concs):
            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

            if contributes_to_signal is not None:
                contributes = contributes_to_signal
            else:
                if init_conc == 0:
                    contributes = False
                else:
                    contributes = True

            for measurement in well.measurements:
                measurement.add_to_blank_states(
                    species_id=species.id,
                    contributes_to_signal=contributes,
                )

        print(
            f"Assigned {species.name} with concentrations of"
            f" {init_concs} {conc_unit} to row {row_id}."
        )

    def _handle_wavelength(self) -> float:
        if len(self.measured_wavelengths) == 1:
            return self.measured_wavelengths[0]

        else:
            raise AttributeError(
                "Argument 'wavelength' must be provided. Measured wavelengths are:"
                f" {self.measured_wavelengths}"
            )

    def get_wells(
        self,
        ids: List[str] = None,
        rows: List[str] = None,
        columns: List[int] = None,
        wavelength: int = None,
    ) -> List[Well]:
        if not wavelength:
            wavelength = self._handle_wavelength()

    def _get_well_by_id(self, _id: str) -> List[Well]:
        for well in self.wells:
            if well.id == _id:
                return well

        raise ValueError(f"No well found with id {_id}")

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

    def _get_wells_by_row_id(self, row_id: str) -> List[Well]:
        return [well for well in self.wells if row_id in well.id and well.measurements]

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
        contributes_to_signal: Optional[bool] = None,
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

            if contributes_to_signal is not None:
                contributes = contributes_to_signal
            else:
                if init_conc == 0:
                    contributes = False
                else:
                    contributes = True

            for measurement in well.measurements:
                measurement.add_to_blank_states(
                    species_id=species.id,
                    contributes_to_signal=contributes,
                )

    def assign_species_from_spreadsheet(
        self,
        species: AbstractSpecies,
        conc_unit: str,
        path: str,
        sheet_name: str,
        header: int = 1,
        index_col: int = 0,
    ):
        df = pd.read_excel(
            path, sheet_name=sheet_name, header=header, index_col=index_col
        )

        for well in self.wells:
            row = well.id[0]
            column = int(well.id[1:])
            init_conc = df.loc[row, column]

            if np.isnan(init_conc):
                continue

            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

            if init_conc == 0:
                contributes_to_signal = False
            else:
                contributes_to_signal = True

            for measurement in well.measurements:
                measurement.add_to_blank_states(
                    species_id=species.id,
                    contributes_to_signal=contributes_to_signal,
                )

    def get_well(self, _id: str) -> Well:
        for well in self.wells:
            if well.id == _id:
                return well

        raise ValueError(f"Well {_id} does not exist.")

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
        reactant_standard: Standard = None,
        sort_measurements_by: AbstractSpecies = None,
        ignore_blank_status: bool = False,
        wavelength: int = None,
        path: str = None,
    ) -> "EnzymeML":
        from MTPHandler.ioutils.enzymeml import create_enzymeml

        return create_enzymeml(
            name=name,
            plate=self,
            detected_reactant=detected_reactant,
            reactant_standard=reactant_standard,
            sort_measurements_by=sort_measurements_by,
            ignore_blank_status=ignore_blank_status,
            wavelength=wavelength,
            path=path,
        )

    def blank_species(
        self,
        species: AbstractSpecies,
        wavelength: int,
    ):
        wells = []
        blanking_wells = []
        for well in self.wells:
            if not well._contains_species(species.id):
                continue

            if not wavelength in [
                measurement.wavelength for measurement in well.measurements
            ]:
                continue

            measurement = well.get_measurement(wavelength)

            if measurement.species_contibutes(species.id):
                wells.append(well)

            if measurement.is_blanked_for(species.id):
                blanking_wells.append(well)

        if len(wells) == 0:
            print(f"{species.name} at {wavelength} nm does not contribute to signal.")
            return

        if len(blanking_wells) == 0:
            raise ValueError(
                f"No wells found to calculate absorption contribution of {species.name} at {wavelength} nm."
            )

        # get mapping of concentration to blank wells
        conc_blank_mapping = self._get_conc_blank_mapping(
            wells=blanking_wells, species=species, wavelength=wavelength
        )

        # apply to wells where species is present in respective concentration
        blanked_wells = []
        for well in wells:
            init_conc = well._get_species_condition(species.id).init_conc
            measurement = well.get_measurement(wavelength)
            blank_state = measurement.get_blank_state(species.id)

            if init_conc not in conc_blank_mapping.keys():
                print(
                    f"Well {well.id} was not blanked for initial"
                    f" {species.name} concentration"
                    f" {init_conc}"
                )
                continue

            if not blank_state.contributes_to_signal:
                continue

            measurement.absorptions = [
                absorption - conc_blank_mapping[init_conc]
                for absorption in measurement.absorptions
            ]

            blank_state.contributes_to_signal = False
            blanked_wells.append(well)

        print(f"Blanked {len(blanked_wells)} wells containing {species.name}.")

    def visualize(self, zoom: bool = False, wavelengths: float = None):
        if zoom:
            shared_yaxes = False
        else:
            shared_yaxes = True

        if not isinstance(wavelengths, list):
            wavelengths = [wavelengths]

        fig = make_subplots(
            rows=self.n_rows,
            cols=self.n_columns,
            shared_xaxes=True,
            subplot_titles=self._generate_possible_well_ids(),
            shared_yaxes=shared_yaxes,
        )
        colors = px.colors.qualitative.Plotly

        for well in self.wells:
            for measurement, color in zip(well.measurements, colors):
                if measurement.wavelength not in wavelengths:
                    continue

                fig.add_trace(
                    go.Scatter(
                        x=self.times,
                        y=measurement.absorptions,
                        name=f"{measurement.wavelength} nm",
                        mode="lines",
                        showlegend=False,
                        line=dict(color=color),
                        hovertemplate="%{y:.2f}<br>",
                    ),
                    col=well.x_position + 1,
                    row=well.y_position + 1,
                )

        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

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

    def _generate_possible_well_ids(self) -> List[str]:
        characters = "ABCDEFGH"
        integers = range(1, 13)  # 1 to 12

        sub_char = characters[: self.n_rows]
        sub_int = integers[: self.n_columns]

        # Generate combinations of characters and integers
        combinations = [
            "".join(item) for item in it.product(sub_char, map(str, sub_int))
        ]

        return combinations

    def _get_conc_blank_mapping(
        self, wells: List[Well], species: AbstractSpecies, wavelength: float
    ) -> Dict[float, List[Well]]:
        blank_measurement_mapping = defaultdict(list)
        for well in wells:
            condition = well._get_species_condition(species.id)

            blank_measurement_mapping[condition.init_conc].append(
                well.get_measurement(wavelength).absorptions
            )

        conc_mean_blank_mapping = {}
        for conc, absorptions in blank_measurement_mapping.items():
            mean_absorption = np.nanmean(absorptions)
            print(
                f"Mean absorption of {species.name} at {conc} {condition.conc_unit}: {mean_absorption}"
            )
            conc_mean_blank_mapping[conc] = mean_absorption

        return conc_mean_blank_mapping

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
