import re
import os
from datetime import datetime
from copy import deepcopy
import numpy as np
from collections import defaultdict
from typing import Dict, Generator, List

import pandas as pd


def read_multiskan(
    cls: "Plate",
    path: str,
    time: float | np.ndarray,
    time_unit: str,
    ph: float = None,
    temperature: float = None,
    temperature_unit: str = None,
    volume: float = None,
    volume_unit: str = None,
) -> "Plate":
    # Get date of measurement
    created = datetime.fromtimestamp(os.path.getctime(path)).strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    # Extract temperature from path
    if not temperature:
        TEMP_PATTERN = r"\d{1,3}deg"
        temperature = re.findall(TEMP_PATTERN, path)[0]
        temperature = re.split("(\d+)", temperature)[1]
        if not temperature:
            raise ValueError("Could not find pH in path. Please specify 'ph'.")

    if not temperature_unit:
        temperature_unit = "C"

    # Extract pH from path
    if not ph:
        PH_PATTERN = r"pH\d+\.\d+"
        ph = re.findall(PH_PATTERN, path)[0]
        ph = float(re.search(r"\d+\.\d+", ph).group())
        if not ph:
            raise ValueError("Could not find pH in path. Please specify pH.")

    # Read file
    data_dict = extract_data(pd.read_csv(path).reset_index())

    # Extract plate dimensions and number of measured timepoints
    n_rows, n_columns, n_timepoints = next(iter(data_dict.values())).shape

    if isinstance(time, list):
        time = np.array([time])

    if len(time) != n_timepoints:
        raise ValueError(
            f"Number of timepoints in data set ({n_timepoints}) does not match "
            f"number of timpoints in provided 'time' array ({len(time)})."
        )

    # Create plate
    plate = cls(
        temperatures=[temperature],
        temperature_unit=temperature_unit,
        times=time,
        time_unit=time_unit,
        date_measured=created,
        n_rows=n_rows,
        n_columns=n_columns,
        measured_wavelengths=list(data_dict.keys()),
    )

    # Add wells to plate
    for wavelength, data in data_dict.items():
        for row_id, row in enumerate(data):
            for column_id, column in enumerate(row):
                id = _coordinates_to_id(column_id, row_id)
                if id not in [well.id for well in plate.wells]:
                    plate.add_to_wells(
                        id=id,
                        x_position=column_id,
                        y_position=row_id,
                        ph=ph,
                    )

                well = [well for well in plate.wells if well.id == id][0]
                well.add_to_measurements(
                    wavelength=wavelength,
                    wavelength_unit="nm",
                    absorptions=column.tolist(),
                )

    return plate


def extract_data(df: pd.DataFrame) -> Dict[int, List[List[float]]]:
    # Get slices of the data corresponding each iteration of measurement
    wavelength_df_dict = _get_plate_dfs(df)

    # Extract data from each iteration
    wavelength_data_dict = defaultdict(list)
    for wavelength, dfs in wavelength_df_dict.items():
        for df in dfs:
            data_of_iteration = df.apply(_extract_row_data, axis=1)
            data_of_iteration = np.array(data_of_iteration.to_list())

            wavelength_data_dict[wavelength].append(data_of_iteration)

    # Convert data to numpy arrays
    for wavelength, data in wavelength_data_dict.items():
        wavelength_data_dict[wavelength] = np.array(data).swapaxes(0, 2).swapaxes(0, 1)

    return wavelength_data_dict


def _get_plate_dfs(df: pd.DataFrame) -> Dict[int, List[pd.DataFrame]]:
    wavelength_slices_dict = defaultdict(list)

    # Get data by wavelength
    wavelengths = df.apply(_get_wavelengths, axis=1).dropna().astype(int).tolist()

    wavelength_dfs = _segment_dataframe(df, "Wavelength")

    for wavelength, wavelength_df in zip(wavelengths, wavelength_dfs):
        iteration_dfs = _segment_dataframe(wavelength_df, "Iteration")

        for iteration_df in iteration_dfs:
            wavelength_slices_dict[wavelength].append(iteration_df)

    return wavelength_slices_dict


def _extract_row_data(row: pd.Series):
    value = row.values[1]
    return [float(x) for x in value.split("\t") if len(x) != 0]


def _segment_dataframe(df: pd.DataFrame, key: str):
    def filter_function(row: pd.Series):
        row_id, value = row.values
        if key in value:
            return row_id

    slice_ids = df.apply(filter_function, axis=1).dropna().astype(int).values.tolist()

    list_of_slices = _list_to_slices(slice_ids)

    return (df.loc[slc] for slc in list_of_slices)


def _list_to_slices(_list: List[int]) -> List[slice]:
    slices = []
    for idx, element in enumerate(_list):
        try:
            start_data_slice = element + 1
            end_data_slice = _list[idx + 1] - 1
            slices.append(slice(start_data_slice, end_data_slice))

        except IndexError:
            start_data_slice = element + 1
            end_data_slice = None
            slices.append(slice(start_data_slice, end_data_slice))

            return slices


def _get_wavelengths(row: pd.Series):
    value = row.values[1]

    if "Wavelength" in value:
        return re.findall(r"Wavelength:(.*)", value)[0]


def _coordinates_to_id(x: int, y: int) -> str:
    return f"{chr(y + 65)}{x+1}"


def id_to_xy(well_id: str):
    return ord(well_id[0].upper()) - 65, int(well_id[1:]) - 1
