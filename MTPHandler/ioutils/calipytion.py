from typing import List
import numpy as np

from CaliPytion.core import Standard, Sample, SignalType
from CaliPytion.tools import Calibrator
from MTPHandler.core.abstractspecies import AbstractSpecies
from MTPHandler.core.well import Well
from MTPHandler.core.protein import Protein


def _get_standard_wells(
    plate: "Plate",
    species: AbstractSpecies,
    wavelength: int,
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

        measurement = well.get_measurement(wavelength)
        if measurement.is_blanked_for(species.id):
            standard_wells.append(well)

        # Add wells with zero concentration to standard wells
        if all(
            [
                blank_state.contributes_to_signal == False
                for blank_state in measurement.blank_states
            ]
        ):
            standard_wells.append(well)

    return standard_wells


def map_to_standard(
    plate: "Plate",
    species: AbstractSpecies,
    wavelength: int = None,
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
        condition = well._get_species_condition(species.id)
        measurement = well.get_measurement(wavelength)
        samples.append(
            Sample(
                id=well.id,
                concentration=condition.init_conc,
                conc_unit=condition.conc_unit,
                signal=float(np.nanmean(measurement.absorptions)),
            )
        )
        phs.append(well.ph)

    # Check if all samples have the same pH
    if not all([ph == phs[0] for ph in phs]):
        raise ValueError(
            f"Samples of standard {species.name} have different pH values: {phs}"
        )
    ph = phs[0]

    # Create standard
    return Standard(
        species_id=species.id,
        name=species.name,
        wavelength=wavelength,
        signal_type=signal_type,
        samples=samples,
        ph=ph,
        temperature=plate.temperatures[0],
        temperature_unit=plate.temperature_unit,
        created=plate.date_measured,
    )


def initialize_calibrator(
    plate: "Plate",
    species: AbstractSpecies,
    wavelength: int = None,
    signal_type: SignalType = SignalType.ABSORBANCE,
    cutoff: float = None,
) -> Calibrator:
    standard = map_to_standard(
        plate=plate,
        species=species,
        wavelength=wavelength,
        signal_type=signal_type,
    )

    return Calibrator.from_standard(standard, cutoff=cutoff)
