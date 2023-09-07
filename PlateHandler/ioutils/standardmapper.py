from typing import List
from collections import defaultdict
from CaliPytion.core.standard import Standard
from PlateHandler.core.species import Species

from PlateHandler.core.well import Well


def wells_to_standard(
        wells: List[Well],
        species: Species
) -> Standard:

    # check if all wells contain the species
    if not all(well._contains_species(species) for well in wells):
        raise AttributeError(
            f"All wells must contain species {species.name}."
        )

    # identify replicates
    conc_well_id_dict = defaultdict(list)
    for well in wells:
        condition = well._get_species_condition(species)
        conc_well_id_dict[condition.init_conc].append(well.id)

    return conc_well_id_dict
