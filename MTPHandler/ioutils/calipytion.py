from __future__ import annotations

from typing import List, Optional

import numpy as np
from calipytion.model import Sample, SignalType, Standard, UnitDefinition
from calipytion.tools.calibrator import Calibrator

from MTPHandler.dataclasses import Molecule, Plate, Protein, Species, Well
from MTPHandler.tools import (
    get_measurement,
    get_species_condition,
    measurement_is_blanked_for,
    well_contains_species,
)


def _get_standard_wells(
    plate: Plate,
    species: Species | Molecule | Protein,
    wavelength: float,
) -> List[Well]:
    # Subset of wells, that contain specified species, do not contain a protein, and are blanked
    protein_ids = [
        species.id for species in plate.species if isinstance(species, Protein)
    ]

    standard_wells = []
    for well in plate.wells:
        if not well_contains_species(well, species.id):
            continue

        if any(
            [well_contains_species(well, catalyst_id) for catalyst_id in protein_ids]
        ):
            continue

        measurement = get_measurement(well, wavelength)
        if measurement_is_blanked_for(measurement, species.id):
            standard_wells.append(well)

        # Add wells with zero concentration to standard wells
        if all(
            [
                blank_state.contributes_to_signal is False
                for blank_state in measurement.blank_states
            ]
        ):
            standard_wells.append(well)

    return standard_wells


def map_to_standard(
    plate: Plate,
    species: Species | Molecule | Protein,
    wavelength: float,
    signal_type: SignalType = SignalType.ABSORBANCE,
) -> Standard:
    standard_wells = _get_standard_wells(
        plate=plate,
        species=species,
        wavelength=wavelength,
    )

    # Map wells to samples of a standard
    samples = []
    phs = []
    for well in standard_wells:
        condition = get_species_condition(well, species.id)
        measurement = get_measurement(well, wavelength)

        samples.append(
            Sample(
                # id=well.id,
                concentration=condition.init_conc,
                conc_unit=UnitDefinition(**condition.conc_unit.model_dump()),
                signal=float(np.nanmean(measurement.absorption)),
            )
        )
        phs.append(well.ph)

    # Check if all samples have the same pH
    if not all([ph == phs[0] for ph in phs]):
        raise ValueError(
            f"Samples of standard {species.name} have different pH values: {phs}"
        )
    ph = phs[0]

    if species.ld_id:
        species_id = species.ld_id
    else:
        species_id = species.id

    temp_unit = UnitDefinition(**plate.temperature_unit.model_dump())

    # Create standard
    return Standard(
        molecule_id=species_id,
        molecule_symbol=species.id,
        molecule_name=species.name,
        wavelength=wavelength,
        signal_type=signal_type,
        samples=samples,
        ph=ph,
        temperature=plate.temperatures[0],
        temp_unit=temp_unit,
    )


def initialize_calibrator(
    plate: Plate,
    wavelength: float,
    species: Species | Molecule | Protein,
    signal_type: SignalType = SignalType.ABSORBANCE,
    cutoff: Optional[float] = None,
) -> Calibrator:
    """
    Initialize a calibrator for a given species.

    Parameters
    ----------
    plate : Plate
        The plate containing the wells.
    wavelength : float
        The wavelength of the measurements.
    species : Species | Molecule | Protein
        The species to calibrate.
    signal_type : SignalType, optional
        The type of signal to calibrate, by default SignalType.ABSORBANCE.

    Returns
    -------
    Calibrator
        The initialized calibrator for creating a calibration curve.

    """
    standard = map_to_standard(
        plate=plate,
        species=species,
        wavelength=wavelength,
        signal_type=signal_type,
    )

    return Calibrator.from_standard(standard, cutoff=cutoff)
