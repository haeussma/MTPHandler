from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional

# from calipytion.core import Standard
import pyenzyme as pe
from pyenzyme.model import DataTypes

from MTPHandler.dataclasses import Molecule, Plate, Protein, Species, Well
from calipytion.model import Standard

import numpy as np

from MTPHandler.tools import get_measurement


class Plate_to_EnzymeMLDocument:

    def __init__(
        self,
        plate: Plate,
        measured_molecule: Molecule,
        standard: Optional[Standard] = None,
    ) -> None:

        assert measured_molecule or standard, (
            "Either a measured molecule or a standard must be provided."
        )

        self.plate = plate
        self.measured_molecule = measured_molecule
        self.standard = standard


    @property
    def temperature(self) -> float:
        return np.mean(self.plate.temperatures).tolist()


    @staticmethod
    def map_protein(protein: Protein) -> pe.Protein:
        return pe.Protein(**protein.model_dump())

    @staticmethod
    def map_small_molecule(molecule: Molecule) -> pe.SmallMolecule:
        return pe.SmallMolecule(
            id=molecule.id,
            ld_id=molecule.ld_id,
            name=molecule.name,
            canonical_smiles=molecule.smiles,
        )

    def map_measurement_data(self, well: Well, wavelength: float) -> list[pe.MeasurementData]:

        measurement = pe.Measurement(
            id=well.id,
            name=f"{self.measured_molecule.name}", #type: ignore
            ph=well.ph,
            temperature=self.temperature,
        )

        photo_meas = get_measurement(well, wavelength)
        data_type = pe.DataTypes.ABSORPTION if not self.standard else pe.DataTypes.CONCENTRATION

        for meas in well.init_conditions:
            if meas.species_id == self.measured_molecule.id:
                is_measured_species = True
            else:
                is_measured_species = False

            measurement.add_to_species(
                species_id=meas.species_id,
                init_conc=meas.init_conc,
                data_type=DataTypes.ABSORPTION if is_measured_species else DataTypes.CONCENTRATION,
                data_unit=photo_meas.data_unit,
                data=photo_meas.absorption if is_measured_species else [],
                time=photo_meas.time if is_measured_species else [],
                time_unit=photo_meas.time_unit,
            )


def enzememl_from_plate(
    name: str,
    plate: Plate,
    detected_molecule: Molecule | Species,
    wavelength: float,
    standard: Standard,
    sort_measurements_by: Optional[Molecule | Protein] = None,
    ignore_blank_status: bool = False,
    path: Optional[str] = None,
) -> pe.EnzymeMLDocument:
    """
    Creates an EnzymeML document based on the information of a `Plate`.


    Args:
        name (str): Name of the `EnzymeMLdocument`
        plate (Plate): Plate with defined and blanked species
        reactant (Reactant): Reactant which was observed during experiment
        protein (Protein): Protein catalyst
        wavelength (int, optional): Detection wavelength of the experiment.
        Defaults to None.

    Raises:
        ValueError: If multiple wavelengths were measured, the wavelength must be specified.

    Returns:
        EnzymeML.EnzymeMLDocument: `EnzymeMLDocument` with defined 'measurements',
        'reactants', 'proteins', and a 'vessel'.
    """

    # Get defined reactants and proteins
    molecules = [
        species for species in plate.species if isinstance(species, (Molecule, Species))
    ]
    proteins = [species for species in plate.species if isinstance(species, Protein)]

    small_mol_dict = {}
    for mol_id, molecule in enumerate(molecules):
        species_id = f"s{molecule.id}"
        small_mol_dict[species_id] = pe.SmallMolecule(
            ld_id=molecule.ld_id,
            name=molecule.name,
            id=species_id,
        )

    # Create measurements
    measurements = assamble_measurements(
        plate=plate,
        wavelength=wavelength,
        detected_reactant=detected_molecule,
        reactant_standard=standard,
        ignore_blank_status=ignore_blank_status,
        proteins=proteins,
        time_unit=plate.time_unit,
        time=plate.times,
    )

    enzymeml = pe.EnzymeMLDocument(
        name=name,
        created=datetime.now().isoformat(),
        measurements=measurements,
        proteins=proteins,
    )

    if sort_measurements_by:
        enzymeml = sort_measurements(enzymeml, sort_measurements_by)

    if path:
        write_document(enzymeml, path)

    return enzymeml


def assamble_measurements(
    plate: "Plate",
    wavelength: float,
    detected_reactant: Molecule | Species,
    reactant_standard: Standard,
    proteins: List[Protein],
    ignore_blank_status: bool,
    time_unit: pe.UnitDefinition,
    time: List[float],
) -> list[pe.Measurement]:
    """
    Creates a list of measurements based on the information of a `Plate`.
    Species are added to the respective measurements based on their initial concentration.
    Data is added to the respective species, which was observed during the experiment.

    Args:
        plate (Plate): `Plate with experimental data`
        wavelength (int): Detection wavelength for the reactant
        reactant (Reactant): Observed `Reactant` species of the experiment
        protein (Protein): Protein catalyst of the reaction

    Raises:
        ValueError: If not all wells with initial conditions have equal time arrays
        ValueError: If not all wells with initial conditions have equal time units
        ValueError: If not all wells with initial conditions have equal concentration units

    Returns:
        List[EnzymeML.Measurement]: List of `Measurement` objects
    """

    # Get wells that contain the specified reactant, contain a catalyst and are blanked
    reaction_wells = get_catalyzed_wells(
        plate=plate,
        wavelength=wavelength,
        detected_reactant=detected_reactant,
        proteins=proteins,
        ignore_blank_status=ignore_blank_status,
    )

    # Group wells by their initial conditions
    well_groups = group_wells_by_init_conditions(wells=reaction_wells)

    # Create measurements
    measurements = []
    for meas_id, well_group in enumerate(well_groups):
        measurement = pe.Measurement(
            id=f"m{meas_id}",
            name=f"{detected_reactant.name} measurement",
            ph=well_group[0].ph,
            temperature=plate.temperatures[0],
            temperature_unit=plate.temperature_unit,
        )

        measurement.species = get_measurement_species(
            measurement=measurement,
            wells=well_group,
            detected_reactant=detected_reactant,
            reactant_standard=reactant_standard,
            time_unit=time_unit,
            time=time,
            wavelength=wavelength,
            ignore_blank_status=ignore_blank_status,
        )

        measurements.append(measurement)

    # Add data to replicates

    return measurements


def group_wells_by_init_conditions(wells: List[Well]) -> List[List[Well]]:
    """Groups wells by their initial conditions."""

    well_mapping = defaultdict(list)
    for well in wells:
        condition_set = frozenset(
            {
                condition.species_id: condition.init_conc
                for condition in well.init_conditions
            }.items()
        )

        well_mapping[condition_set].append(well)

    return [wells for wells in well_mapping.values()]


def get_measurement_species(
    measurement: EnzymeML.Measurement,
    wells: List[Well],
    detected_reactant: Reactant,
    reactant_standard: Standard,
    time_unit: str,
    time: List[float],
    wavelength: float,
    ignore_blank_status: bool,
) -> List[EnzymeML.MeasurementData]:
    """
    Creates a list of `MeasurementData` objects for a `Measurement` object.

    Args:
        measurement (EnzymeML.Measurement): `Measurement` object
        wells (List[Well]): List of `Well` objects, that form a replicate
        reactant (Reactant): Observed `Reactant` species of the experiment

    Returns:
        List[EnzymeML.MeasurementData]: List of `MeasurementData` objects
    """

    measurement_datas = []
    for condition in wells[0].init_conditions:
        measurement_data = EnzymeML.MeasurementData(
            init_conc=condition.init_conc,
            unit=condition.conc_unit,
            measurement_id=measurement.id,
            species_id=condition.species_id,
        )

        if condition.species_id == detected_reactant.id:
            measurement_data.replicates = get_replicates(
                measurement_data=measurement_data,
                wells=wells,
                standard=reactant_standard,
                time_unit=time_unit,
                time=time,
                wavelength=wavelength,
                ignore_blank_status=ignore_blank_status,
            )

        measurement_datas.append(measurement_data)

    return measurement_datas


def get_replicates(
    measurement_data: EnzymeML.MeasurementData,
    wells: List[Well],
    standard: Standard,
    time_unit: str,
    time: List[float],
    wavelength: float,
    ignore_blank_status: bool,
) -> List[EnzymeML.Replicate]:
    """
    Creates a list of `Replicate` objects for a `MeasurementData` object.

    Args:
        measurement_data (EnzymeML.MeasurementData): `MeasurementData` object
        wells (List[Well]): List of `Well` objects, that form a replicate

    Returns:
        List[EnzymeML.Replicate]: List of `Replicate` objects
    """

    replicates = []
    if standard:
        units = [sample.conc_unit for sample in standard.samples]
        if not all([u == units[0] for u in units]):
            raise ValueError("Concentration units of standard samples are not equal.")
        unit = str(units[0])

        for well in wells:
            data = well.get_measurement(wavelength)
            replicate = EnzymeML.Replicate(
                id=well.id,
                species_id=measurement_data.species_id,
                measurement_id=measurement_data.measurement_id,
                data_type=DataTypes.CONCENTRATION,
                data_unit=unit,
                time_unit=time_unit,
                time=time,
                data=data.to_concentration(
                    standard=standard, ignore_blank_status=ignore_blank_status
                ),
            )
            replicates.append(replicate)

    else:
        for well in wells:
            replicate = EnzymeML.Replicate(
                id=well.id,
                species_id=measurement_data.species_id,
                measurement_id=measurement_data.measurement_id,
                data_type=DataTypes.ABSORPTION,
                data_unit="dimensionless",
                time_unit=time_unit,
                time=time,
                data=well.get_measurement(wavelength).absorptions,
            )
            replicates.append(replicate)

    return replicates


def get_replicates_mapping(
    wells: List[Well], mapping_reactant: Reactant
) -> Dict[float, List[str]]:
    """
    Creates a mapping of initial concentrations to well ids.
    Keys are unique initial concentrations of the observed `Reactant`.
    Values are lists of well ids, forming a `Replicate`.

    Args:
        wells (List[Well]): List of `Well` objects, that form a replicate
        reactant (Reactant): Observed `Reactant` species of the experiment

    Returns:
        Dict[float, List[str]]: Mapping of initial concentrations to well ids
    """

    replicates_mapping = defaultdict(list)
    for well in wells:
        condition = well._get_species_condition(mapping_reactant.id)
        replicates_mapping[condition.init_conc].append(well)

    return replicates_mapping


def get_catalyzed_wells(
    plate: "Plate",
    wavelength: int,
    detected_reactant: Reactant,
    proteins: List[Protein],
    ignore_blank_status: bool,
) -> List[Well]:
    """
    Returns a list of wells, that contain the specified species,
    contain a catalyst, and are blanked.

    Args:
        plate (Plate): `Plate` object with experimental data
        wavelength (int): Detection wavelength for the reactant
        reactant (Reactant): Observed `Reactant` species of the experiment
        protein (Protein): Protein catalyst of the reaction

    Returns:
        List[Well]: List of `Well` objects
    """

    # Subset of wells, that contain one or more proteins, and one or more reactants that are `constant=False`

    catalyzed_wells = []
    for well in plate.wells:
        if not any([well._contains_species(protein.id) for protein in proteins]):
            continue

        if not ignore_blank_status:
            if not well.get_measurement(wavelength).is_blanked_for(
                detected_reactant.id
            ):
                continue
        catalyzed_wells.append(well)

    print(f"Found {len(catalyzed_wells)} catalyzed wells")

    return catalyzed_wells


def sort_measurements(
    enzymeMLDocument: EnzymeML.EnzymeMLDocument,
    sort_species: AbstractSpecies,
) -> EnzymeML.EnzymeMLDocument:
    measurements_dict = {}
    for measurement in enzymeMLDocument.measurements:
        for species in measurement.species:
            if species.species_id == sort_species.id:
                measurements_dict[species.init_conc] = measurement

    sorted_keys = list(sorted(measurements_dict.keys()))
    sorted_measurements = [measurements_dict[key] for key in sorted_keys]

    enzymeMLDocument.measurements = sorted_measurements

    return enzymeMLDocument


def write_document(enzymeMLDocument: EnzymeML.EnzymeMLDocument, path: str):
    """Writes file to specified path.
    Supported formats are `json` and `yaml`.

    Args:
        enzymeMLDocument (EnzymeML.EnzymeMLDocument): `EnzymeMLDocument` object
        path (str): Path to file

    Raises:
        ValueError: If path is not ending with `json` or `yaml`
    """

    format = path.split(".")[-1]

    if format == "json":
        outfile = enzymeMLDocument.json()
    elif format == "yaml":
        outfile = enzymeMLDocument.yaml()
    else:
        raise ValueError(f"Format {format} is not supported.")

    with open(path, "w") as f:
        f.write(outfile)
