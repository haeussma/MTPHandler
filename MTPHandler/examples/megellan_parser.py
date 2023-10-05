from email.policy import default
import re
import os
from datetime import datetime
from copy import deepcopy
import numpy as np
from collections import defaultdict
from typing import Dict, Generator, List

import pandas as pd


def parse_magellan(
    cls: "Plate",
    path: str,
    ph: float,
    wavelengths: List[float],
):
    created = datetime.fromtimestamp(os.path.getctime(path)).strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    df = pd.read_excel(path, header=None)

    # Define the format of the input datetime string
    date_format = "%A, %B %d, %Y: %H:%M:%S"

    data = defaultdict(list)
    temperatures = []
    times = []
    dates = []
    for row in df.iterrows():
        timecourser_data = row[1].values[0]
        if not isinstance(timecourser_data, str):
            break
        else:
            date_str, time_str, temperature_str = timecourser_data.split("/")
            temp_value, temp_unit = temperature_str.strip().split("°")
            temperatures.append(float(temp_value))
            time, time_unit = time_str[1:-1].split(" ")

            times.append(time)
            time_unit = time_unit
            dates.append(datetime.strptime(date_str.strip(), date_format))

    df = df.dropna()

    for row in df.iterrows():
        for element in row[1].values:
            if isinstance(element, str):
                key = element
            else:
                data[key].append(element)

    n_rows = 8
    n_columns = 12

    plate = cls(
        ph=ph,
        created=created,
        n_rows=n_rows,
        n_columns=n_columns,
        measured_wavelengths=wavelengths,
        temperature_unit=temp_unit,
    )

    for well_id, abso_list in data.items():
        x_pos, y_pos = id_to_xy(well_id)
        plate.add_to_wells(
            id=well_id,
            absorption=abso_list,
            time=times,
            x_position=x_pos,
            y_position=y_pos,
            wavelength=wavelengths[0],
        )

    return plate


def id_to_xy(well_id: str):
    return int(well_id[1:]) - 1, ord(well_id[0].upper()) - 65
