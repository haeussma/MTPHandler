import time
from typing import List
from sdRDM import DataModel

from PlateHandler.modified.plate import Plate
from PlateHandler.modified.species import Species

EnzymeML = DataModel.from_markdown(
    "test_scheme/enzymeml.md")

# TODO: from urllib.request import request try: url open


def create_measurements(
        plate: Plate,
        reactant: Species,
        catalyst: Species,
        wavelength: int = None
) -> List["EnzymeML.Measurement"]:

    # Get wells that contain the specified reactant, contain a catalyst and are blanked
    reaction_wells = get_catalyzed_wells(
        plate=plate,
        wavelength=wavelength,
        reactant=reactant,
        catalyst=catalyst,
    )

    # Check consistency of time arrays and time units
    if not all([well.time_unit == reaction_wells[0].time_unit for well in reaction_wells]):
        raise ValueError("Time units of wells are not equal.")
    time_unit = reaction_wells[0].time_unit

    if not all([well.time == reaction_wells[0].time for well in reaction_wells]):
        raise ValueError("Time values of wells are not equal.")
    time = reaction_wells[0].time

    # Create measurements
    measurements = []
    for well in reaction_wells:
        measurement = EnzymeML.Measurement(
            name=well.id,
            global_time=time,
            global_time_unit=time_unit,
            ph=plate.ph,
            temperature=plate.temperature,
            temperature_unit=plate.temperature_unit,
        )
        measurements.append(measurement)

    return measurements


def get_catalyzed_wells(
        plate: Plate,
        wavelength: int,
        reactant: Species,
        catalyst: Species,
) -> List["EnzymeML.Measurement"]:

    # handel wavelength
    if not wavelength:
        if len(plate.measured_wavelengths) == 1:
            wavelength = plate.measured_wavelengths[0]
        else:
            raise AttributeError(
                f"Argument 'wavelength' must be provided. Measured wavelengths are: {plate.measured_wavelengths}"
            )

    # Subset of wells, that contain specified species, contain a catalyst and are blanked
    reaction_wells = []
    for well in plate.wells:
        if not well._contains_species(reactant):
            continue

        if not well._contains_species(catalyst):
            continue

        if not well._is_blanked(reactant):
            continue

        reaction_wells.append(well)

    return reaction_wells
