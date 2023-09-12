from typing import List
from collections import defaultdict

from CaliPytion.tools.standardcurve import StandardCurve
from CaliPytion.core.analyte import Analyte
from CaliPytion.core.standard import Standard
from CaliPytion.core.series import Series
import numpy as np
from MTPHandler.modified.species import Species
from MTPHandler.modified.plate import Plate
from MTPHandler.modified.well import Well


def _get_standard_wells(
        plate: Plate,
        species: Species,
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

    # Subset of wells, that contain specified species, do not contain a catalyst and are blanked
    standard_wells = []
    for well in wells:
        if not well._contains_species(species):
            continue

        if plate._get_catalyst():
            catalyst = plate._get_catalyst()
            if well._contains_species(catalyst):
                continue

        if well._is_blanked(species):
            standard_wells.append(well)

    return standard_wells


def create_standard_curve(
        plate: Plate,
        species: Species,
        wavelength: int = None,
        cutoff_signal: float = None
) -> Analyte:

    standard_wells = _get_standard_wells(
        plate=plate,
        species=species,
        wavelength=wavelength,
    )

    # group wells by same initial concentration
    conc_well_id_map = defaultdict(list)
    units = []
    for well in standard_wells:
        condition = well._get_species_condition(species)
        conc_well_id_map[condition.init_conc].append(well.id)
        units.append(condition.conc_unit)

    if all(unit == units[0] for unit in units):
        unit = units[0]

    concentrations = list(conc_well_id_map.keys())
    absorption_replicates = list(map(list, zip(*conc_well_id_map.values()))
                                 )  # transpose list of well ids

    # create series
    series = []
    for replicate in absorption_replicates:
        wells = [plate.get_well(well_id) for well_id in replicate]
        series.append(
            Series(values=[np.mean(well.absorption) for well in wells]))

    # create Standard
    standard = Standard(
        wavelength=wavelength,
        concentration=concentrations,
        absorption=series,
        concentration_unit=unit
    )

    # create Analyte
    analyte = Analyte(
        name=species.name,
        ph=plate.ph,
        temperature=plate.temperature,
        temperature_unit=plate.temperature_unit,
        standard=[standard],
    )

    return StandardCurve.from_analyte(
        analyte=analyte,
        wavelength=wavelength,
        cutoff_signal=cutoff_signal,
        blank_data=False,
    )
