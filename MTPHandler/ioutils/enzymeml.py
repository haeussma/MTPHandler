from typing import Dict, List
from sdRDM import DataModel
from collections import defaultdict
from datetime import datetime

from CaliPytion import Standard

from MTPHandler.modified.initcondition import InitCondition
from MTPHandler.modified.reactant import Reactant
from MTPHandler.modified.protein import Protein
from MTPHandler.modified.well import Well

# Specify EnzymeML version
URL = "https://github.com/EnzymeML/enzymeml-specifications.git"
COMMIT = "72c3d8be4a094983667a7aa62fb599fbc9f7351c"

EnzymeML = DataModel.from_git(URL, COMMIT)
SBOTerm = EnzymeML.enums.SBOTerm
DataTypes = EnzymeML.enums.DataTypes


def create_enzymeml(
        name: str,
        plate: "Plate",
        detected_reactant: Reactant,
        reactant_standard: Standard,
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
    reactants = [
        species for species in plate.species if isinstance(species, Reactant)]

    proteins = [species for species in plate.species if isinstance(
        species, Protein)]

    # Create measurements
    measurements = create_measurements(
        plate=plate,
        wavelength=wavelength,
        detected_reactant=detected_reactant,
        reactant_standard=reactant_standard,
        proteins=proteins
    )

    # Create EnzymeMLDocument
    enzymeMLDocument = EnzymeML.EnzymeMLDocument(
        name=name,
        created=plate.created,
        vessels=[plate._define_dummy_vessel()],
        reactants=reactants,
        proteins=proteins,
        measurements=measurements
    )

    if path:
        write_doument(enzymeMLDocument, path)

    return enzymeMLDocument


def write_doument(
        enzymeMLDocument: EnzymeML.EnzymeMLDocument,
        path: str
):
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


def create_measurements(
        plate: "Plate",
        wavelength: int,
        detected_reactant: Reactant,
        reactant_standard: Standard,
        proteins: List[Protein]
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
        proteins=proteins
    )

    # Check consistency of time arrays and time units
    if not all([well.time_unit == reaction_wells[0].time_unit for well in reaction_wells]):
        raise ValueError("Time units of wells are not equal.")
    time_unit = reaction_wells[0].time_unit

    if not all([well.time == reaction_wells[0].time for well in reaction_wells]):
        raise ValueError("Time values of wells are not equal.")
    time = reaction_wells[0].time

    if not all(well._get_species_condition(detected_reactant.id).conc_unit == reaction_wells[0]._get_species_condition(detected_reactant.id).conc_unit for well in reaction_wells):
        raise ValueError("Concentration units of wells are not equal.")
    conc_unit = reaction_wells[0]._get_species_condition(
        detected_reactant.id).conc_unit

    # Get mapping of replicates
    replicates_mapping = get_replicates_mapping(
        wells=reaction_wells,
        detected_reactant=detected_reactant,
    )

    # Create measurements
    measurements = []
    for init_conc, well_ids in replicates_mapping.items():
        measurement = EnzymeML.Measurement(
            name=f"{detected_reactant.name} {init_conc} {conc_unit}",
            global_time=time,
            global_time_unit=time_unit,
            ph=plate.ph,
            temperature=plate.temperature,
            temperature_unit=plate.temperature_unit,
        )

        replicate_wells = [plate.get_well(well_id) for well_id in well_ids]
        measurement.species = get_measurement_species(
            measurement=measurement,
            wells=replicate_wells,
            detected_reactant=detected_reactant,
            reactant_standard=reactant_standard,
        )

        measurements.append(measurement)

    # Add data to replicates

    return measurements


def get_measurement_species(
        measurement: EnzymeML.Measurement,
        wells: List[Well],
        detected_reactant: Reactant,
        reactant_standard: Standard
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

    # TODO: Assumes all wells have equal list of `InitCondition` objects. Might not be the case.

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
                standard=reactant_standard
            )

        measurement_datas.append(measurement_data)

    return measurement_datas


def get_replicates(
        measurement_data: EnzymeML.MeasurementData,
        wells: List[Well],
        standard: Standard
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
            raise ValueError(
                "Concentration units of standard samples are not equal.")
        unit = str(units[0])

        for well in wells:
            replicate = EnzymeML.Replicate(
                id=well.id,
                species_id=measurement_data.species_id,
                measurement_id=measurement_data.measurement_id,
                data_type=DataTypes.CONCENTRATION,
                data_unit=unit,
                time_unit=well.time_unit,
                time=well.time,
                data=well.to_concentration(standard=standard),
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
                time_unit=well.time_unit,
                time=well.time,
                data=well.absorption,
            )
            replicates.append(replicate)

    return replicates


def get_replicates_mapping(
        wells: List[Well],
        detected_reactant: Reactant
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
        condition = well._get_species_condition(detected_reactant.id)
        replicates_mapping[condition.init_conc].append(well.id)

    return replicates_mapping


def get_catalyzed_wells(
        plate: "Plate",
        wavelength: int,
        detected_reactant: Reactant,
        proteins: List[Protein]
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
    for well in plate.get_wells(wavelength=wavelength):

        if not any([well._contains_species(protein.id) for protein in proteins]):
            continue

        if not well._contains_species(detected_reactant.id):
            continue

        if not well._is_blanked_for(detected_reactant.id):
            continue

        catalyzed_wells.append(well)

    print(f"Found {len(catalyzed_wells)} catalyzed wells")

    return catalyzed_wells
