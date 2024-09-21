import importlib.resources as pkg_resources

import httpx
import toml

from mtphandler.model import InitCondition, PhotometricMeasurement, Plate, Well


def read_static_file(path, filename: str):
    """Reads a static file from the specified library path.

    Args:
        path (Module): Import path of the library.
        filename (str): The name of the file to read.

    Returns:
        dict: The contents of the file as a dictionary.
    """

    source = pkg_resources.files(path).joinpath(filename)
    with pkg_resources.as_file(source) as file:
        return toml.load(file)


def get_measurement(well: Well, wavelength: float) -> PhotometricMeasurement:
    """
    Get the measurement object for a given well and wavelength.

    Args:
        well (Well): The well object.
        wavelength (float): The wavelength of the measurement.

    Returns:
        PhotometricMeasurement: The measurement object.

    Raises:
        ValueError: If no measurement is found for the given well and wavelength.
    """

    for measurement in well.measurements:
        if measurement.wavelength == wavelength:
            return measurement

    raise ValueError(
        f"No measurement found for well {well.id} at wavelength {wavelength}."
    )


def well_contains_species(
    well: Well, species_id: str, conc_above_zero: bool = False
) -> bool:
    """Check if a well contains a species with the given ID, and optionally, if its concentration is above zero.

    Args:
        well (Well): The well to check.
        species_id (str): The ID of the species.
        conc_above_zero (bool): If True, checks if the species' concentration is above zero.

    Returns:
        bool: True if the species is present in the well (and has a concentration above zero if conc_above_zero is True), otherwise False.
    """
    for condition in well.init_conditions:
        if condition.species_id == species_id:
            # If conc_above_zero is True, check if concentration is > 0
            if conc_above_zero:
                return condition.init_conc > 0
            # Otherwise, just return True if species is present
            return True

    return False


def handle_blank_status(
    well: Well,
    species_id: str,
    init_conc: float,
    contributes_to_signal: bool | None,
):
    """Add blank status to the measurements of a well.
    If the concentration is 0, the species does not contribute to the signal.
    If the concentration is not 0, the species contributes to the signal unless
    overwriten by the `contributes_to_signal` argument.

    Args:
        well (Well): Well for which to add blank status.
        species_id (str): ID of the species.
        init_conc (float): Initial concentration of the species.
        contributes_to_signal (bool | None): Whether the species contributes to the signal.
    """
    if contributes_to_signal is None:
        if init_conc == 0:
            contributes = False
        else:
            contributes = True
    else:
        contributes = contributes_to_signal

    for measurement in well.measurements:
        measurement.add_to_blank_states(
            species_id=species_id,
            contributes_to_signal=contributes,
        )


def measurement_is_blanked_for(
    measurement: PhotometricMeasurement, target_id: str
) -> bool:
    """Checks if a the measurement is blanked for a given species target."""

    target_contributes = None
    others_contribute = []

    for state in measurement.blank_states:
        if state.species_id == target_id:
            target_contributes = state.contributes_to_signal

        else:
            others_contribute.append(state.contributes_to_signal)

    if target_contributes is None:
        raise ValueError(f"Species {target_id} not found in blank states")

    return target_contributes and not any(others_contribute)


def get_species_condition(well: Well, species_id: str) -> InitCondition:
    for condition in well.init_conditions:
        if condition.species_id == species_id:
            return condition

    raise ValueError(f"Species {species_id} not found in well {well.id}")


def pubchem_request_molecule_name(pubchem_cid: int) -> str:
    """Retrieves molecule name from PubChem database based on CID."""

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/Title/JSON"
    response = httpx.get(url)

    if response.status_code == 200:
        res_dict = response.json()
        try:
            molecule_name = res_dict["PropertyTable"]["Properties"][0]["Title"]
            return molecule_name
        except (KeyError, IndexError):
            raise ValueError(
                "Unexpected response structure while retrieving molecule name from PubChem"
            )
    else:
        raise ValueError("Failed to retrieve molecule name from PubChem")


def get_well(plate: Plate, well_id: str) -> Well:
    for well in plate.wells:
        if well.id.lower() == well_id.lower():
            return well

    raise ValueError(f"Well {well_id} not found")
