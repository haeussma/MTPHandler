from __future__ import annotations

import itertools as it
import re
import warnings
from collections import defaultdict
from datetime import datetime as Datetime
from types import NoneType
from typing import Dict, List, Literal, Optional, Union
from uuid import uuid4

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sdRDM
from calipytion.core import SignalType
from calipytion.tools.calibrator import Calibrator
from lxml.etree import _Element
from plotly.subplots import make_subplots
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from rich.jupyter import print
from sdRDM.base.datatypes import Identifier, Unit
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from MTPHandler.ioutils import initialize_calibrator
from MTPHandler.readers.factory import MTPReaderFactory

from .initcondition import InitCondition
from .photometricmeasurement import PhotometricMeasurement
from .protein import Protein
from .reactant import Reactant
from .species import Species
from .well import Well


class Plate(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Description of a microtiter plate consisting of wells."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    n_rows: int = element(
        description="Number of rows on the plate",
        tag="n_rows",
        json_schema_extra=dict(),
    )

    n_cols: int = element(
        description="Number of columns on the plate",
        tag="n_cols",
        json_schema_extra=dict(),
    )

    name: Optional[str] = element(
        description="Name of the plate",
        default=None,
        tag="name",
        json_schema_extra=dict(),
    )

    date_measured: Optional[Datetime] = element(
        description="Date and time when the plate was measured",
        default=None,
        tag="date_measured",
        json_schema_extra=dict(),
    )

    times: List[float] = element(
        description=(
            "Time points of the measurement, corresponding to temperature measurements"
        ),
        default_factory=ListPlus,
        tag="times",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    time_unit: Optional[Unit] = element(
        description="Unit of the time",
        default=None,
        tag="time_unit",
        json_schema_extra=dict(),
    )

    temperatures: List[float] = element(
        description="Thermostat temperature",
        default_factory=ListPlus,
        tag="temperatures",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    temperature_unit: Optional[Unit] = element(
        description="Unit of the temperature",
        default=None,
        tag="temperature_unit",
        json_schema_extra=dict(),
    )

    wells: List[Well] = element(
        description="List of wells on the plate",
        default_factory=ListPlus,
        tag="wells",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    species: List[Species] = element(
        description="List of species present in wells of the plate",
        default_factory=ListPlus,
        tag="species",
        json_schema_extra=dict(
            multiple=True,
        ),
    )

    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/MTPHandler"
    )
    _commit: Optional[str] = PrivateAttr(
        default="c1ace6fef38751e79952a37557221f0540ca9d77"
    )

    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                isinstance(i, _Element) for i in value
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)

        return self

    def add_to_wells(
        self,
        x_pos: int,
        y_pos: int,
        ph: Optional[float] = None,
        init_conditions: List[InitCondition] = ListPlus(),
        measurements: List[PhotometricMeasurement] = ListPlus(),
        volume: Optional[float] = None,
        volume_unit: Optional[Unit] = None,
        id: Optional[str] = None,
        **kwargs,
    ) -> Well:
        """
        This method adds an object of type 'Well' to attribute wells

        Args:
            id (str): Unique identifier of the 'Well' object. Defaults to 'None'.
            x_pos (): X position of the well on the plate.
            y_pos (): Y position of the well on the plate.
            ph (): pH of the reaction. Defaults to None
            init_conditions (): List of initial conditions of different species. Defaults to ListPlus()
            measurements (): List of photometric measurements. Defaults to ListPlus()
            volume (): Volume of the reaction. Defaults to None
            volume_unit (): Unit of the volume. Defaults to None
        """

        params = {
            "x_pos": x_pos,
            "y_pos": y_pos,
            "ph": ph,
            "init_conditions": init_conditions,
            "measurements": measurements,
            "volume": volume,
            "volume_unit": volume_unit,
        }

        if id is not None:
            params["id"] = id

        obj = Well(**params)

        self.wells.append(obj)

        return self.wells[-1]

    def add_to_species(
        self,
        name: Optional[str] = None,
        references: List[Identifier] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> Species:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            name (): Name of the species. Defaults to None
            references (): List of references to the species. Defaults to ListPlus()
        """

        params = {
            "name": name,
            "references": references,
        }

        if id is not None:
            params["id"] = id

        obj = Species(**params)

        self.species.append(obj)

        return self.species[-1]

    def add_reactant_to_species(
        self,
        name: Optional[str] = None,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        references: List[Identifier] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> Reactant:
        """
        This method adds an object of type 'Reactant' to attribute species

        Args:
            id (str): Unique identifier of the 'Reactant' object. Defaults to 'None'.
            name (): Name of the species. Defaults to None
            smiles (): SMILES representation of the species. Defaults to None
            inchi (): InChI representation of the species. Defaults to None
            references (): List of references to the Reactant. Defaults to ListPlus()
        """

        params = {
            "name": name,
            "smiles": smiles,
            "inchi": inchi,
            "references": references,
        }

        if id is not None:
            params["id"] = id

        obj = Reactant(**params)

        self.species.append(obj)

        return self.species[-1]

    def add_protein_to_species(
        self,
        name: Optional[str] = None,
        sequence: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[Identifier] = None,
        references: List[Identifier] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> Protein:
        """
        This method adds an object of type 'Protein' to attribute species

        Args:
            id (str): Unique identifier of the 'Protein' object. Defaults to 'None'.
            name (): Name of the species. Defaults to None
            sequence (): Amino acid sequence of the protein. Defaults to None
            organism (): Organism the protein originates from. Defaults to None
            organism_tax_id (): NCBI taxonomy ID of the organism. Defaults to None
            references (): List of references to the protein. Defaults to ListPlus()
        """

        params = {
            "name": name,
            "sequence": sequence,
            "organism": organism,
            "organism_tax_id": organism_tax_id,
            "references": references,
        }

        if id is not None:
            params["id"] = id

        obj = Protein(**params)

        self.species.append(obj)

        return self.species[-1]

    def add_protein(
        self,
        name: Optional[str] = None,
        sequence: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[Identifier] = None,
        references: List[Identifier] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> Protein:
        params = {
            "id": id,
            "name": name,
            "sequence": sequence,
            "organism": organism,
            "organism_tax_id": organism_tax_id,
            "references": references,
            **kwargs,
        }

        if id in [species.id for species in self.species]:
            old_species = self.get_species(id)
            del self.species[self.species.index(old_species)]
        protein = self.add_protein_to_species(**params)
        protein.add_object_term("https://schema.org/Protein")

        return protein

    def add_species(
        self,
        name: Optional[str] = None,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        references: List[Identifier] = ListPlus(),
        id: Optional[str] = None,
        **kwargs,
    ) -> Reactant:
        """Add a reactant to the species list.

        Args:
            name (Optional[str], optional): Name of the reactant. Defaults to None.
            smiles (Optional[str], optional): Smiles code of the molecule. Defaults to None.
            inchi (Optional[str], optional): Inchi code of the molecule. Defaults to None.
            references (List[Identifier], optional): An link pointing to a reference. Defaults to ListPlus().
            id (Optional[str], optional): ID of the reactant. This can be e.g. an Inchi-Key
                or Chebi identifier. Defaults to None.

        Returns:
            Reactant: The reactant object.
        """

        if not id:
            warnings.warn(
                """No ID provided for species. This might cause issues when referencing the species.
                Consider providing an ID which references the species such as a Chebi ID."""
            )

        params = {
            "id": id,
            "name": name,
            "smiles": smiles,
            "inchi": inchi,
            "references": references,
            **kwargs,
        }

        if id in [species.id for species in self.species]:
            old_species = self.get_species(id)
            del self.species[self.species.index(old_species)]
        reactant = self.add_to_species(**params)

        return reactant

    def assign_species(
        self,
        species: Species,
        init_conc: Union[float, List[float]],
        conc_unit: str,
        to: Literal["all", "rows", "columns", "all except"],
        ids: Union[str, List[str], int, List[int]] = None,
        contributes_to_signal: bool | None = None,
    ):
        """
        Assigns a species to specific wells on the plate based on the provided criteria.

        Args:
            species (Species): The species object to assign to the wells.
            init_conc (Union[float, List[float]]): The initial concentration(s) of the species.
            conc_unit (str): The unit of concentration.
            to (Literal["all", "rows", "columns", "all except"]): The target location(s) for assigning the species.
                This can be individual or multiple wells as well as individual or multiple rows or columns.
            ids (Union[str, List[str], int, List[int]], optional): The specific row or column IDs where the species should be assigned.
                Defaults to None.
            contributes_to_signal (bool, optional): Indicates if the assigned species contributes to the signal.
                Defaults to None.

        Raises:
            AttributeError: If the 'species' argument does not reference an `id` of a Species listed in `species`.
            AttributeError: If the 'to' argument is not one of ["all", "rows", "columns", "all_except"].

        Returns:
            None
        """
        # Handle species
        if isinstance(species, str):
            species = self.get_species(species)
        elif isinstance(species, Species):
            pass
        else:
            raise AttributeError(
                """Argument 'species' must reference an `id` of a Species listed in `species`."""
            )

        cases = ["rows", "columns", "all", "all except"]
        if to not in cases:
            raise AttributeError(f"Argument 'to' must be one of {cases}.")

        if not isinstance(init_conc, list):
            init_conc = [init_conc]

        if not isinstance(ids, list) and not isinstance(ids, NoneType):
            ids = [ids]

        if to == "all":
            if isinstance(init_conc, list):
                if len(init_conc) == 1:
                    init_conc = init_conc[0]

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
            if isinstance(init_conc, list):
                if len(init_conc) == 1:
                    init_conc = init_conc[0]

            self.assign_species_to_all_except(
                well_ids=ids,
                species=species,
                init_conc=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

    def assign_species_to_all(
        self,
        species: Species,
        init_conc: float,
        conc_unit: Unit,
        contributes_to_signal: bool | None,
    ):
        for well in self.wells:
            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

            self._handle_blank_status(
                well, species.id, init_conc, contributes_to_signal
            )

        print(
            f"Assigned [bold magenta]{species.name}[/] ({species.id}) with"
            f" {init_conc} {conc_unit} to all wells."
        )

    def assign_species_to_columns(
        self,
        column_ids: list[int],
        species: Species,
        init_concs: list[float],
        conc_unit: str,
        contributes_to_signal: bool | None,
    ):
        # Handle column_ids
        if not all([isinstance(column_id, int) for column_id in column_ids]):
            raise AttributeError("Argument 'column_ids' must be a list of integers.")

        if not all([column_id <= self.n_cols for column_id in column_ids]):
            raise AttributeError(
                "Argument 'column_ids' must be a list of integers between 1 and"
                f" {self.n_cols+1}."
            )

        columns = []
        for column_id in column_ids:
            wells = [well for well in self.wells if well.x_pos + 1 == column_id]
            wells = sorted(wells, key=lambda x: x.y_pos)
            columns.append(wells)

        # Handle init_concs
        if len(init_concs) == 1:
            init_concs = init_concs * len(columns[0])

        for wells in columns:
            assert len(init_concs) == len(wells), f"""
            Number of initial concentrations ({len(init_concs)}) does not match number
            of wells ({len(wells)}) in columns ({column_ids}).
            """

            for well, init_conc in zip(wells, init_concs):
                well.add_to_init_conditions(
                    species_id=species.id,
                    init_conc=init_conc,
                    conc_unit=conc_unit,
                )

                self._handle_blank_status(
                    well, species.id, init_conc, contributes_to_signal
                )

        print(
            f"Assigned [bold magenta]{species.name}[/] ({species.id}) with"
            f" concentrations of {init_concs} {conc_unit} to columns {column_ids}."
        )

    def assign_species_to_rows(
        self,
        row_ids: List[str],
        species: Species,
        init_concs: List[float],
        conc_unit: str,
        contributes_to_signal: bool | None,
    ):
        # Handle row_ids

        if isinstance(row_ids, str):
            row_ids = [row_ids]

        if not all([isinstance(row_id, str) for row_id in row_ids]):
            raise AttributeError("Argument 'row_ids' must be a list of strings.")

        rows = []
        for row_id in row_ids:
            wells = [well for well in self.wells if row_id in well.id]
            wells = sorted(wells, key=lambda x: x.x_pos)
            rows.append(wells)

        for wells in rows:
            assert len(init_concs) == len(wells), f"""
            Number of initial concentrations ({len(init_concs)}) does not match number
            of wells ({len(wells)}) in rows ({row_ids}).
            """

            for well, init_conc in zip(wells, init_concs):
                well.add_to_init_conditions(
                    species_id=species.id,
                    init_conc=init_conc,
                    conc_unit=conc_unit,
                )

                self._handle_blank_status(
                    well, species.id, init_conc, contributes_to_signal
                )

        print(
            f"Assigned [bold magenta]{species.name}[/] ({species.id}) with"
            f" {init_conc} {conc_unit} to rows {row_ids}."
        )

    @staticmethod
    def _handle_blank_status(
        well: Well,
        species_id: str,
        init_conc: float,
        contributes_to_signal: bool | None,
    ):
        if contributes_to_signal is None:
            if init_conc == 0:
                contributes = False
            else:
                contributes = True
        else:
            contributes = contributes_to_signal

        for measurement in well.measurements:
            measurement.add_to_blank_states(
                species_id=species_id,
                contributes_to_signal=contributes,
            )

    def assign_species_to_all_except(
        self,
        well_ids: List[str],
        species: Species,
        init_conc: float,
        conc_unit: str,
        contributes_to_signal: bool | None,
    ):
        # validate well_id
        self._validate_well_id(well_ids)

        wells = (well for well in self.wells if well.id not in well_ids)
        for well in wells:
            well.add_to_init_conditions(
                species_id=species.id,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

            self._handle_blank_status(
                well, species.id, init_conc, contributes_to_signal
            )

        print(
            f"Assigned [bold magenta]{species.name}[/] ({species.id}) with"
            f" {init_conc} {conc_unit} to all wells except {well_ids}."
        )

    def assign_species_from_spreadsheet(
        self,
        species: Species,
        conc_unit: str,
        path: str,
        sheet_name: str | None = None,
        header: int = 0,
        index: int = 0,
        contributes_to_signal: bool | None = None,
    ):
        """Allows to read in a spreadsheet harboring a plate map describing the
        initial concentrations of a species in the wells of the plate.
        The parser expects the first row to contain the row ids numbered from 1 to 12.
        The first column should contain the column ids labeled from A to H.
        If additional columns or rows are present, adjust the `header` and `index`
        arguments accordingly.

        Args:
            species (Species): Species object to assign initial concentrations to.
            conc_unit (str): Unit of the concentration.
            path (str): Path to the spreadsheet.
            sheet_name (str | None, optional): Name of the sheet containing the plate map.
                Defaults to None.
            header (int, optional): Number of columns before the header starts.
                Defaults to 0.
            index (int, optional): Number of columns before the index column starts.
                Defaults to 0.
            contributes_to_signal (bool | None, optional): Whether the species contributes
                to the signal. Defaults to None.
        """
        df = pd.read_excel(path, sheet_name=sheet_name, header=header, index_col=index)

        if isinstance(species, str):
            species = self.get_species(species)

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

            self._handle_blank_status(
                well, species.id, init_conc, contributes_to_signal
            )

        print(
            f"Assigned initial concentrations for [bold magenta]{species.name}[/]"
            f" ({species.id}) from {path}."
        )

    def get_well(self, _id: str) -> Well:
        for well in self.wells:
            if well.id == _id:
                return well

        raise ValueError(f"Well {_id} does not exist.")

    def calibrate(
        self,
        species: Union[Reactant, Protein],
        wavelength: int,
        signal_type: SignalType,
        cutoff: float = None,
    ) -> Calibrator:
        return initialize_calibrator(
            plate=self,
            species=species,
            wavelength=wavelength,
            signal_type=signal_type,
            cutoff=cutoff,
        )

    # def to_enzymeml(
    #     self,
    #     name: str,
    #     detected_reactant: Reactant,
    #     reactant_standard: Standard = None,
    #     sort_measurements_by: AbstractSpecies = None,
    #     ignore_blank_status: bool = False,
    #     wavelength: int = None,
    #     path: str = None,
    # ) -> "EnzymeML":
    #     return create_enzymeml(
    #         name=name,
    #         plate=self,
    #         detected_reactant=detected_reactant,
    #         reactant_standard=reactant_standard,
    #         sort_measurements_by=sort_measurements_by,
    #         ignore_blank_status=ignore_blank_status,
    #         wavelength=wavelength,
    #         path=path,
    #     )

    def blank_species(
        self,
        species: Species,
        wavelength: int,
    ):
        wells = []
        blanking_wells = []
        for well in self.wells:
            if not well._contains_species(species.id):
                continue

            if wavelength not in [
                measurement.wavelength for measurement in well.measurements
            ]:
                continue

            measurement = well.get_measurement(wavelength)

            if measurement.species_contibutes(species.id):
                wells.append(well)

            if measurement.is_blanked_for(species.id):
                blanking_wells.append(well)

        if len(wells) == 0:
            print(
                f"{species.name} ({species.id}) at {wavelength} nm does not contribute"
                " to signal. Is the specified wavelength correct?print"
            )

        if len(blanking_wells) == 0:
            raise ValueError(
                "No wells found to calculate absorption contribution of"
                f" {species.name} ({species.id}) at {wavelength} nm."
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

            measurement.absorption = [
                absorption - conc_blank_mapping[init_conc]
                for absorption in measurement.absorption
            ]

            blank_state.contributes_to_signal = False
            blanked_wells.append(well)

        print(f"Blanked {len(blanked_wells)} wells containing {species.name}.\n")

    def visualize(
        self, zoom: bool = False, wavelengths: list[float] = [], static: bool = False
    ):
        if zoom:
            shared_yaxes = False
        else:
            shared_yaxes = True

        if not wavelengths:
            wavelengths = [self.wells[0].measurements[0].wavelength]

        if not isinstance(wavelengths, list):
            wavelengths = [wavelengths]

        fig = make_subplots(
            rows=self.n_rows,
            cols=self.n_cols,
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
                        y=measurement.absorption,
                        name=f"{measurement.wavelength} nm",
                        mode="lines",
                        showlegend=False,
                        line=dict(color=color),
                        hovertemplate="%{y:.2f}<br>",
                    ),
                    col=well.x_pos + 1,
                    row=well.y_pos + 1,
                )

        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        fig.update_layout(
            plot_bgcolor="white",
            hovermode="x",
            title=dict(
                text=f"{self.temperatures[0]} {self.temperature_unit.name}",
            ),
            margin=dict(l=20, r=20, t=100, b=20),
        )

        if static:
            fig.show("svg")

        fig.show()

    def get_species(self, _id: str) -> Union[Reactant, Protein]:
        for species in self.species:
            if species.id == _id:
                return species

        raise ValueError(f"No species found with id {_id}")

    def _generate_possible_well_ids(self) -> List[str]:
        characters = "ABCDEFGH"
        integers = range(1, 13)  # 1 to 12

        sub_char = characters[: self.n_rows]
        sub_int = integers[: self.n_cols]

        # Generate combinations of characters and integers
        combinations = [
            "".join(item) for item in it.product(sub_char, map(str, sub_int))
        ]

        return combinations

    def _get_conc_blank_mapping(
        self, wells: List[Well], species: Species, wavelength: float
    ) -> Dict[float, List[Well]]:
        blank_measurement_mapping = defaultdict(list)
        for well in wells:
            condition = well._get_species_condition(species.id)

            blank_measurement_mapping[condition.init_conc].append(
                well.get_measurement(wavelength).absorption
            )

        conc_mean_blank_mapping = {}
        for conc, absorptions in blank_measurement_mapping.items():
            mean_absorption = np.nanmean(absorptions)
            std_absorption = np.nanstd(absorptions)
            std_perc = abs(std_absorption / mean_absorption) * 100
            print(
                f"Mean absorption of [bold magenta]{species.name}[/] ({species.id}) at"
                f" {conc} {condition.conc_unit.name}: {mean_absorption:.4f} ±"
                f" {std_perc:.0f}%  calculated based on wells"
                f" {[well.id for well in wells]}."
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
    def from_file(
        cls: Plate,
        path: str,
        ph: Optional[float] = None,
        wavelength: Optional[float] = None,
        time: Optional[List[float]] = None,
        time_unit: Optional[str] = None,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
    ) -> Plate:
        return MTPReaderFactory.read(
            cls=cls,
            path=path,
            ph=ph,
            wavelength=wavelength,
            time=time,
            time_unit=time_unit,
            temperature=temperature,
            temperature_unit=temperature_unit,
        )

    @classmethod
    def from_reader(cls):
        warnings.warn(
            "Method 'from_reader' is deprecated. Use 'from_file' instead.",
            DeprecationWarning,
        )
