from __future__ import annotations

from datetime import datetime

import pandas as pd
from typing_extensions import Optional

from MTPHandler.model import Plate
from MTPHandler.readers.utils import id_to_xy
from MTPHandler.units import C, second


def read_tekan_spark(
    path: str,
    ph: Optional[float],
) -> Plate:
    df = pd.read_excel(path)

    cycle_no_row_index = df[df.iloc[:, 0].str.contains("Cycle Nr.", na=False)].index[0]
    meta_df = df.iloc[:cycle_no_row_index, :]
    data_df = df.iloc[cycle_no_row_index:, :].reset_index(drop=True)

    meta_df = (
        meta_df.dropna(how="all")
        .dropna(axis=1, how="all")
        .set_index(meta_df.columns[0])
    )
    time_measured = meta_df.loc["Start Time"].dropna(axis=1, how="all").values[0][0]
    time_measured = datetime.strptime(time_measured, "%Y-%m-%d %H:%M:%S")

    wavelength = meta_df.loc["Measurement wavelength"].dropna().iloc[0]

    data_df = data_df.set_index(data_df.columns[0])
    column_names = data_df.iloc[0, :].tolist()
    data_df.columns = column_names
    data_df = data_df[1:].reset_index(drop=True).dropna(axis=1, how="all")
    first_nan_index = data_df.isna().any(axis=1).idxmax()
    data_df = data_df.iloc[:first_nan_index, :]

    time_series = data_df.pop("Time [s]") / 60
    temp_series = data_df.pop("Temp. [°C]")

    plate = Plate(
        date_measured=str(time_measured),
        temperatures=temp_series.values.tolist(),
        temperature_unit=C,
        time_unit=second,
        times=time_series.values.tolist(),
    )

    for column in data_df.columns:
        x, y = id_to_xy(column)
        well = plate.add_to_wells(
            id=column,
            x_pos=x,
            y_pos=y,
            ph=ph,
        )

        well.add_to_measurements(
            wavelength=wavelength,
            wavelength_unit="nm",
            absorption=data_df[column].values.tolist(),
            time=time_series.values.tolist(),
            time_unit=second,
        )

    return plate


if __name__ == "__main__":
    from devtools import pprint

    path = "/Users/max/Documents/GitHub/MTPHandler/docs/examples/data/tekan_spark.xlsx"
    from MTPHandler.model import Plate

    p = read_tekan_spark(path, 7.4)

    pprint(p.wells[0])
