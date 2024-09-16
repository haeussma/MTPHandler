import importlib.resources as pkg_resources

import toml

from MTPHandler.dataclasses import InitCondition, PhotometricMeasurement, Well


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


def well_contains_species(well: Well, species_id: str) -> bool:
    """Check if a well contains a species with the given id."""
    for condition in well.init_conditions:
        if condition.species_id == species_id:
            return True

    return False


def handle_blank_status(
    well: Well,
    species_id: str,
    init_conc: float,
    contributes_to_signal: bool | None,
):
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
    measurement: PhotometricMeasurement, species_id: str
) -> bool:
    """Checks if the measurement is blanked for a given species."""

    if species_id not in [state.species_id for state in measurement.blank_states]:
        raise ValueError(f"Species {species_id} is not present in this well.")

    target_contributes = [
        state.contributes_to_signal
        for state in measurement.blank_states
        if state.species_id == species_id
    ][0]

    others_contribute = [
        state.contributes_to_signal
        for state in measurement.blank_states
        if state.species_id != species_id
    ]

    if target_contributes and not any(others_contribute):
        return True

    return False


def get_species_condition(well: Well, species_id: str) -> InitCondition:
    for condition in well.init_conditions:
        if condition.species_id == species_id:
            return condition

    raise ValueError(f"Species {species_id} not found in well {well.id}")
