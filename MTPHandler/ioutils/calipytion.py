from os import name
from typing import List
from collections import defaultdict
import numpy as np

from CaliPytion import Calibrator, Standard
from CaliPytion.modified.sample import Sample
from MTPHandler.core.abstractspecies import AbstractSpecies
from MTPHandler.core.well import Well
from MTPHandler.core.protein import Protein


def _get_standard_wells(
    plate: "Plate", species: AbstractSpecies, wavelength: int
) -> List[Well]:
    # handel wavelength
    if not wavelength:
        if len(plate.measured_wavelengths) == 1:
            wavelength = plate.measured_wavelengths[0]
        else:
            raise AttributeError(
                f"Argument 'wavelength' must be provided. Measured wavelengths are: {plate.measured_wavelengths}"
            )

    # Subset of wells, that contain specified species, do not contain a protein, and are blanked
    protein_ids = [
        species.id for species in plate.species if isinstance(species, Protein)
    ]

    standard_wells = []
    for well in plate.wells:
        if not well._contains_species(species.id):
            continue

        if any([well._contains_species(catalyst_id) for catalyst_id in protein_ids]):
            continue

        if well._get_measurement(wavelength).is_blanked_for(species.id):
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
        condition = well._get_species_condition(species.id)
        measurement = well.get_measurement(wavelength)
        samples.append(
            Sample(
                id=well.id,
                concentration=condition.init_conc,
                conc_unit=condition.conc_unit,
                signal=np.nanmean(measurement.absorptions),
            )
        )

    # Create standard
    return Standard(
        species_id=species.id,
        name=species.name,
        wavelength=wavelength,
        samples=samples,
        ph=plate.ph,
        temperature=plate.temperatures[0],
        temperature_unit=plate.temperature_unit,
        created=plate.date_measured,
    )


def initialize_calibrator(
    plate: "Plate",
    species: AbstractSpecies,
    wavelength: int = None,
    cutoff: float = None,
) -> Calibrator:
    standard = map_to_standard(
        plate=plate,
        species=species,
        wavelength=wavelength,
    )

    return Calibrator.from_standard(standard, cutoff=cutoff)
