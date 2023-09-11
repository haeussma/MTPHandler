from typing import List
from collections import defaultdict


from CaliPytion.core.standard import Standard
from PlateHandler.core.species import Species
from PlateHandler.core.plate import Plate
from PlateHandler.core.well import Well


def create_standard(
        plate: Plate,
        species: Species,
        wavelength: int = None,
) -> Standard:

    # handel wavelength
    if not wavelength:
        if len(plate.measured_wavelengths) == 1:
            wavelength = plate.measured_wavelengths[0]
        else:
            raise AttributeError(
                f"Argument 'wavelength' must be provided. Measured wavelengths are: {plate.measured_wavelengths}"
            )

    wells = plate.get_wells(wavelength=wavelength)
    print(f"all wells: {len(wells)}")

    # get wells, of species if no enzyme is present and 'was_blanked' is True
    # for all species except the one of interest.

    standard_wells = []
    for well in wells:
        if not well._contains_species(species):
            print(f"no species {well.id, species.id}")
            continue

        if plate._get_catalyst_id():
            print("yes has catalyst")
            catalyst_id = plate._get_catalyst_id()
            if well._contains_species(catalyst_id):
                print(f"catalyst {well.id}")
                continue

        if well._contains_catalyst():
            print(f"catalyst {well.id}")
            continue
        if well._all_blanked_except_species(species):
            standard_wells.append(well)

    return standard_wells

    # # check if all wells contain the species
    # wells = [well for well in wells if well.contains_species(species)]

    # # group wells by same initial concentration
    # conc_well_id_dict = defaultdict(list)
    # for well in wells:
    #     condition = well._get_species_condition(species)
    #     conc_well_id_dict[condition.init_conc].append(well.id)

    # concentrations = list(conc_well_id_dict.keys())
    # absorption_replicate_ids = list(map(list, zip(*conc_well_id_dict.values()))
    #                    )  # transpose list of well ids

    # # create standard
    # standard = Standard(
    #     wavelength=
