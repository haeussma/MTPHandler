import itertools as it
import re
from collections import defaultdict
from datetime import datetime as Datetime
from types import NoneType
from typing import Callable, Dict, List, Literal, Optional, Union
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
from sdRDM.base.datatypes import Identifier, Unit
from sdRDM.base.listplus import ListPlus
from sdRDM.tools.utils import elem2dict

from MTPHandler.ioutils import initialize_calibrator

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

        protein = self.add_protein_to_species(**params)
        protein.add_object_term("https://schema.org/Protein")

        return protein

    def add_reactant(
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
        params = {
            "id": id,
            "name": name,
            "smiles": smiles,
            "inchi": inchi,
            "references": references,
            **kwargs,
        }

        reactant = self.add_reactant_to_species(**params)

        return reactant

    def assign_species(
        self,
        species_id: str,
        init_conc: Union[float, List[float]],
        conc_unit: str,
        to: Literal["all", "rows", "columns", "except"],
        ids: Union[str, List[str], int, List[int]] = None,
        contributes_to_signal: Optional[bool] = None,
    ):
        # Handle species
        if isinstance(species_id, str):
            species_id = self.get_species(species_id)
        elif isinstance(species_id, Protein) or isinstance(species_id, Reactant):
            species_id = species_id.id
        else:
            raise AttributeError(
                """Argument 'species_id' must reference an `if` of a species listed in `species` 
                attribute of a `Plate` or be an instance of Protein or Reactant."""
            )

        cases = ["rows", "columns", "all", "except"]
        if to not in cases:
            raise AttributeError(f"Argument 'to' must be one of {cases}.")

        if not isinstance(init_conc, list):
            init_conc = [init_conc]

        if not isinstance(ids, list) and not isinstance(ids, NoneType):
            ids = [ids]

        if to == "all":
            self.assign_species_to_all(
                species_id=species_id,
                init_conc=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        elif to == "columns":
            self.assign_species_to_columns(
                column_ids=ids,
                species_id=species_id,
                init_concs=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        elif to == "rows":
            self.assign_species_to_rows(
                row_ids=ids,
                species_id=species_id,
                init_concs=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        else:
            self.assign_species_to_all_except(
                well_ids=ids,
                species_id=species_id,
                init_conc=init_conc,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

        return

    def assign_species_to_all(
        self,
        species_id: str,
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
                species_id=species_id.id,
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
                    species_id=species_id.id,
                    contributes_to_signal=contributes,
                )

        print(f"Assigned {species_id.name} to all wells.")

    def assign_species_to_columns(
        self,
        column_ids: List[int],
        species_id: str,
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
                        species_id=species_id.id,
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
                            species_id=species_id.id,
                            contributes_to_signal=contributes,
                        )

        print(
            f"Assigned {species_id.name} with concentrations of"
            f" {init_concs} {conc_unit} to columns {column_ids}."
        )

    def assign_species_to_rows(
        self,
        row_ids: List[str],
        species_id: str,
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
                species_id=species_id,
                conc_unit=conc_unit,
                contributes_to_signal=contributes_to_signal,
            )

    def _assign_species_to_row(
        self,
        row_id: str,
        init_concs: List[float],
        species_id: str,
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
                species_id=species_id.id,
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
                    species_id=species_id.id,
                    contributes_to_signal=contributes,
                )

        print(
            f"Assigned {species_id.name} with concentrations of"
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
        species_id: str,
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
                species_id=species_id.id,
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
                    species_id=species_id.id,
                    contributes_to_signal=contributes,
                )

    def assign_species_from_spreadsheet(
        self,
        species_id: str,
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
                species_id=species_id,
                init_conc=init_conc,
                conc_unit=conc_unit,
            )

            if init_conc == 0:
                contributes_to_signal = False
            else:
                contributes_to_signal = True

            for measurement in well.measurements:
                measurement.add_to_blank_states(
                    species_id=species_id,
                    contributes_to_signal=contributes_to_signal,
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
        species_id: str,
        wavelength: int,
    ):
        wells = []
        blanking_wells = []
        for well in self.wells:
            if not well._contains_species(species_id):
                continue

            if wavelength not in [
                measurement.wavelength for measurement in well.measurements
            ]:
                continue

            measurement = well.get_measurement(wavelength)

            if measurement.species_contibutes(species_id):
                wells.append(well)

            if measurement.is_blanked_for(species_id):
                blanking_wells.append(well)

        if len(wells) == 0:
            print(f"{species_id} at {wavelength} nm does not contribute to signal.")
            return

        if len(blanking_wells) == 0:
            raise ValueError(
                "No wells found to calculate absorption contribution of"
                f" {species_id} at {wavelength} nm."
            )

        # get mapping of concentration to blank wells
        conc_blank_mapping = self._get_conc_blank_mapping(
            wells=blanking_wells, species_id=species_id, wavelength=wavelength
        )

        # apply to wells where species is present in respective concentration
        blanked_wells = []
        for well in wells:
            init_conc = well._get_species_condition(species_id.id).init_conc
            measurement = well.get_measurement(wavelength)
            blank_state = measurement.get_blank_state(species_id.id)

            if init_conc not in conc_blank_mapping.keys():
                print(
                    f"Well {well.id} was not blanked for initial"
                    f" {species_id.name} concentration"
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

        print(f"Blanked {len(blanked_wells)} wells containing {species_id.name}.")

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
                text=f"{self.temperatures[0]} °{self.temperature_unit}",
            ),
            margin=dict(l=20, r=20, t=100, b=20),
        )

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
        sub_int = integers[: self.n_columns]

        # Generate combinations of characters and integers
        combinations = [
            "".join(item) for item in it.product(sub_char, map(str, sub_int))
        ]

        return combinations

    def _get_conc_blank_mapping(
        self, wells: List[Well], species_id: str, wavelength: float
    ) -> Dict[float, List[Well]]:
        blank_measurement_mapping = defaultdict(list)
        for well in wells:
            condition = well._get_species_condition(species_id.id)

            blank_measurement_mapping[condition.init_conc].append(
                well.get_measurement(wavelength).absorptions
            )

        conc_mean_blank_mapping = {}
        for conc, absorptions in blank_measurement_mapping.items():
            mean_absorption = np.nanmean(absorptions)
            print(
                f"Mean absorption of {species_id} at {conc} {condition.conc_unit}:"
                f" {mean_absorption}"
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
