from os import name
from typing import List
from collections import defaultdict
import numpy as np

from CaliPytion import Calibrator, Standard
from CaliPytion.modified.sample import Sample
from MTPHandler.modified.abstractspecies import AbstractSpecies
from MTPHandler.modified.well import Well


def _get_standard_wells(
        plate: "Plate",
        species: AbstractSpecies,
        wavelength: int
) -> List[Well]:

    # handel wavelength
    if not wavelength:
        if len(plate.measured_wavelengths) == 1:
            wavelength = plate.measured_wavelengths[0]
        else:
            raise AttributeError(
                f"Argument 'wavelength' must be provided. Measured wavelengths are: {plate.measured_wavelengths}"
            )

    wells = plate.get_wells(wavelength=wavelength)

    # Subset of wells, that contain specified species, do not contain a catalyst, and are blanked
    standard_wells = []
    for well in wells:
        if not well._contains_species(species):
            continue

        if plate._get_catalyst():
            catalyst = plate._get_catalyst()
            if well._contains_species(catalyst):
                continue

        if well._is_blanked(species.id):
            standard_wells.append(well)

    return standard_wells


def map_to_standard(
        plate: "Plate",
        species: AbstractSpecies,
        wavelength: int = None,
) -> Standard:

    standard_wells = _get_standard_wells(
        plate=plate,
        species=species,
        wavelength=wavelength,
    )

    # Map wells to samples of a standard
    samples = []
    for well in standard_wells:
        condition = well._get_species_condition(species)
        samples.append(
            Sample(
                id=well.id,
                concentration=condition.init_conc,
                conc_unit=condition.conc_unit,
                signal=np.mean(well.absorption)
            )
        )

    # Create standard
    return Standard(
        species_id=species.id,
        name=species.name,
        wavelength=wavelength,
        samples=samples,
        ph=plate.ph,
        temperature=plate.temperature,
        temperature_unit=plate.temperature_unit,
        created=plate.created,
    )


def initialize_calibrator(
        plate: "Plate",
        species: AbstractSpecies,
        wavelength: int = None,
        cutoff: float = 3,
) -> Calibrator:

    standard = map_to_standard(
        plate=plate,
        species=species,
        wavelength=wavelength,
    )

    return Calibrator.from_standard(standard, cutoff=cutoff)
