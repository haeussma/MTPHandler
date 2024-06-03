from collections import defaultdict
from datetime import datetime
from typing import Dict, List

from calipytion.core import Standard
from sdRDM import DataModel

from MTPHandler.core.protein import Protein
from MTPHandler.core.reactant import Reactant
from MTPHandler.core.well import Well

# Specify EnzymeML version
URL = "https://github.com/EnzymeML/enzymeml-specifications.git"
COMMIT = "5e5f05b9dc76134305b8f9cef65271e35563ac76"

EnzymeML = DataModel.from_git(URL, COMMIT)
SBOTerm = EnzymeML.enums.SBOTerm
DataTypes = EnzymeML.enums.DataTypes


def create_enzymeml(
    name: str,
    plate: "Plate",
    detected_reactant: Reactant,
    reactant_standard: Standard,
    sort_measurements_by: AbstractSpecies = None,
    ignore_blank_status: bool = False,
    wavelength: int = None,
    path: str = None,
) -> EnzymeML.EnzymeMLDocument:
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

    # handel wavelength
    if not wavelength:
        if len(plate.measured_wavelengths) == 1:
            wavelength = plate.measured_wavelengths[0]
        else:
            raise AttributeError(
                f"Argument 'wavelength' must be provided. Measured wavelengths are: {plate.measured_wavelengths}"
            )

    # Get defined reactants and proteins
    reactants = [species for species in plate.species if isinstance(species, Reactant)]

    proteins = [species for species in plate.species if isinstance(species, Protein)]

    # Create measurements
    measurements = assamble_measurements(
        plate=plate,
        wavelength=wavelength,
        detected_reactant=detected_reactant,
        reactant_standard=reactant_standard,
        ignore_blank_status=ignore_blank_status,
        proteins=proteins,
        time_unit=plate.time_unit,
        time=plate.times,
    )

    enzymeml = EnzymeML.EnzymeMLDocument(
        name=name,
        created=datetime.now(),
        vessels=[plate._define_dummy_vessel()],
        measurements=measurements,
        reactants=reactants,
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
    detected_reactant: Reactant,
    reactant_standard: Standard,
    proteins: List[Protein],
    ignore_blank_status: bool,
    time_unit: str,
    time: List[float],
) -> List[EnzymeML.Measurement]:
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
        measurement = EnzymeML.Measurement(
            id=f"m{meas_id}",
            name=f"{detected_reactant.name} measurement",
            global_time=plate.times,
            global_time_unit=plate.time_unit,
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
