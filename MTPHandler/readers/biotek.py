from __future__ import annotations

import re
from typing import TYPE_CHECKING
from typing_extensions import Optional
import pandas as pd
from datetime import datetime
from datetime import time

from MTPHandler.readers.utils import id_to_xy

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate

def read_biotek(
    cls: Plate,
    path: str,
    ph: float|None = None,
) -> Plate:

    df = pd.read_excel(path)

    date = get_row_by_value(df, "Date")[-1]
    time = get_row_by_value(df, "Time")[-1]
    timestamp = datetime.combine(date, time)

    wavelength_cell = df.iloc[19,1]
    wavelengths = re.findall(r'\d+', wavelength_cell)
    wavelengths = list(map(int, wavelengths))

    row_index_int_map = df.iloc[:, 0].apply(lambda cell: isinstance(cell, int))
    data_block_starts = [index for index, value in enumerate(row_index_int_map) if value]

    for row_index, (block_start, wavelength) in enumerate(zip(data_block_starts, wavelengths)):

        try:
            block = df.iloc[block_start+2: data_block_starts[row_index +1],:]
        except IndexError:
            block = df.iloc[block_start+2:,:]

        block = block.drop("Unnamed: 0", axis=1).reset_index(drop=True)

        all_nan_rows = block.isna().all(axis=1)
        first_all_nan_index = all_nan_rows.idxmax() if all_nan_rows.any() else None

        # drop rows if any of the values are NaN
        block.iloc[:first_all_nan_index, :]
        block = block.dropna(how="any", axis=0)
        column_names = block.iloc[0, :].tolist()
        block.columns = column_names
        block = block[1:].reset_index(drop=True)


        times = block.pop("Time").values
        cleaned_time = ([normalize_datetime_to_reference(t) for t in times])

        # relative time
        start_time = cleaned_time[0]
        relative_time = [float((t - start_time).total_seconds()/60) for t in cleaned_time]

        # Temperature
        temperature = block.pop(column_names[1]).values

        if row_index == 0:
            plate = cls(
                date_measured=timestamp,
                n_rows=8,
                n_cols=12,
                times=relative_time,
                time_unit="min",
                temperatures=temperature,
                temperature_unit="°C",
            )

        for column_name in column_names[2:]:
            x, y = id_to_xy(column_name)

            try:
                well = plate.get_well(column_name)
            except ValueError:
                well = plate.add_to_wells(
                    id=column_name,
                    x_pos=x,
                    y_pos=y,
                    ph=ph
                )

            assert(len(block[column_name].values.tolist()) == len(relative_time)), f"Length of data and time does not match for well {column_name}"

            well.add_to_measurements(
                wavelength=wavelength,
                wavelength_unit="nm",
                absorption=block[column_name].values.tolist(),
                time=relative_time,
            )


    return plate



def get_row_by_value(df: pd.DataFrame, value: str) -> list:
    row_df = df[df.iloc[:, 0].values == value]
    row_df = row_df.reset_index(drop=True)
    row_df = row_df.dropna(axis=1, how="all")
    return row_df.loc[0].values.tolist()

def normalize_datetime_to_reference(dt, reference_date=datetime(1899, 12, 31)):
    # Use the reference date for year and month, keep day, hour, minute, second, and microsecond
    if isinstance(dt, time):
        normalized_dt = datetime(
            reference_date.year,
            reference_date.month,
            reference_date.day,
            dt.hour,
            dt.minute,
            dt.second,
        )
        return normalized_dt

    else:
        return dt




if __name__ == "__main__":
    from MTPHandler.core import Plate

    path = "/Users/max/Documents/GitHub/MTPHandler/tests/data/ BioTek_Epoch2.xlsx"

    print(read_biotek(Plate, path))
