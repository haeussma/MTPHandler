from __future__ import annotations

import logging
import math
import re
from collections import defaultdict
from datetime import datetime
from typing import TYPE_CHECKING, Optional

import pandas as pd

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate

LOGGER = logging.getLogger(__name__)


def read_magellan(
    cls: Plate,
    path: str,
    wavelength: float,
    ph: Optional[float] = None,
):
    df = pd.read_excel(path, header=None)

    # Define the format of the input datetime string
    date_format = "%A, %B %d, %Y: %H:%M:%S"
    WELL_ID_PATTERN = r"[A-H][0-9]{1,2}"

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
            time_unit = time_unit.replace("sec", "s")
            dates.append(datetime.strptime(date_str.strip(), date_format))

    created = dates[0]

    df = df.dropna(how="all")

    for row in df.iterrows():
        first_cell = str(row[1].values[0])
        if not re.findall(WELL_ID_PATTERN, first_cell):
            continue
        for element in row[1].values:
            if isinstance(element, str):
                key = element
            elif math.isnan(element):
                continue
            else:
                data[key].append(element)

    n_rows = 8
    n_columns = 12

    plate = cls(
        date_measured=created,
        n_rows=n_rows,
        n_cols=n_columns,
        temperature_unit=temp_unit,
        temperatures=temperatures,
        time_unit=time_unit,
        times=times,
    )

    for well_id, abso_list in data.items():
        x_pos, y_pos = id_to_xy(well_id)
        well = plate.add_to_wells(
            ph=ph if ph else None,
            id=well_id,
            x_pos=x_pos,
            y_pos=y_pos,
        )
        well.add_to_measurements(
            wavelength=wavelength,
            wavelength_unit="nm",
            absorption=abso_list,
            # blank_states=[],
        )

    return plate


def id_to_xy(well_id: str):
    return int(well_id[1:]) - 1, ord(well_id[0].upper()) - 65
