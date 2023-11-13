import re
from datetime import datetime
from copy import deepcopy
import numpy as np
from typing import List

import pandas as pd
from MTPHandler.core import Well, PhotometricMeasurement


def parse_spectramax(
    cls: "Plate",
    path: str,
    ph: float = None,
):
    df = pd.read_csv(
        path,
        sep="delimiter",
        encoding="utf-16",
        engine="python",
        skiprows=15,
    )

    df = df.map(lambda x: x.split("\t"))

    # Get date of measurement
    last_saved = re.findall(
        r"\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2} [APMapm]{2}", df.iloc[-1, 0][0]
    )[0]
    created = datetime.strptime(last_saved, "%Y/%m/%d %I:%M:%S %p")

    wavelengths = df.iloc[0, 0][-6]
    wavelengths = [
        float(wavelength) for wavelength in wavelengths.split(" ") if wavelength != ""
    ]

    # data
    datas = df.iloc[2:-9]["~End"].tolist()
    time_pattern = r"\d{2}:\d{2}:\d{2}"
    times = []
    temperatures = []
    block_start_ids = []
    for data in datas:
        if re.match(time_pattern, data[0]):
            times.append(to_seconds(data[0]))
            temperatures.append(float(data[1]))
            block_start_ids.append(datas.index(data))

    time_blocks = []
    for index, start_id in enumerate(block_start_ids):
        try:
            time_block = datas[start_id : block_start_ids[index + 1]]
        except IndexError:
            time_block = datas[start_id:]

        time_block[0] = time_block[0][2:]  # remove time and temperature from block
        for row_id, row in enumerate(time_block):
            if "" not in row:
                continue

            wavelength_entry = []
            wavelength_entries = []
            for item in row:
                if item == "":
                    if wavelength_entry:
                        wavelength_entries.append(wavelength_entry)
                        wavelength_entry = []  # reset

                else:
                    wavelength_entry.append(item)

            if item != "":
                wavelength_entries.append(wavelength_entry)

            time_block[row_id] = wavelength_entries

        time_blocks.append(time_block)

    # Swap dimensions: rows, columns, wavelengths, timecourse
    data = np.array(time_blocks).astype(float)
    data = data.swapaxes(0, 3)
    data = data.swapaxes(0, 1)

    # create wells
    wells = []
    for row_id, row in enumerate(data):
        for column_id, column in enumerate(row):
            well = Well(
                id=_coordinates_to_id(column_id, row_id),
                x_position=column_id,
                y_position=row_id,
            )
            for wavelength_id, wavelength in enumerate(column):
                well.add_to_measurements(
                    wavelength=wavelengths[wavelength_id],
                    wavelength_unit="nm",
                    absorptions=wavelength.tolist(),
                )
            wells.append(well)

    # Create plate
    plate = cls(
        n_rows=8,
        n_columns=12,
        date_measured=created,
        times=times,
        time_unit="s",
        temperatures=temperatures,
        temperature_unit="C",
        ph=ph,
        wells=wells,
    )

    return plate

    # # Read raw data, parse metadata

    # df = df.iloc[13:-9].reset_index().drop(columns="index")
    # initial_substrates = [0, 5, 10, 15, 25, 50, 75, 100, 150, 200]
    # time = []
    # array = []

    # # Extract measurement data
    # for index, row in df.iterrows():
    #     if index % 6 == 0:
    #         time.append(row.values[0][0])
    #         df.loc[index, "##BLOCKS= 2"] = row.values[0][2:]

    #     data = row.values[0]
    #     array.append([string_to_float(x) for x in data])

    # # Cleaning and restructuring of data
    # array: np.ndarray = np.array(array)
    # array = array[~np.isnan(array)]
    # array = array.reshape(11, 6, 20)
    # array = array.swapaxes(0, 2)

    # # (wavelength, concentration, control, replicates, data)
    # array = array.reshape((2, 10, 2, 3, 11))

    # # (wavelength, control, concentration, replicates, data)
    # array = array.swapaxes(1, 2)

    # substrate = array[0][0]
    # substrate_control = array[0][1]
    # product = array[1][0]
    # product_control = array[1][1]

    # measurements_product_control = []
    # for data, init_substrate in zip(product_control, initial_substrates):
    #     measurement = {}
    #     measurement["initial_substrate"] = 0
    #     reps = []
    #     for replicate in data:
    #         rep_dict = {}
    #         rep_dict["replicate"] = replicate.tolist()
    #         reps.append(rep_dict)
    #     measurement["data"] = reps

    #     measurements_product_control.append(measurement)

    # measurements_substrate_control = []
    # for data, init_substrate in zip(substrate_control, initial_substrates):
    #     measurement = {}
    #     measurement["initial_substrate"] = init_substrate
    #     reps = []
    #     for replicate in data:
    #         rep_dict = {}
    #         rep_dict["replicate"] = replicate.tolist()
    #         reps.append(rep_dict)
    #     measurement["data"] = reps

    #     measurements_substrate_control.append(measurement)

    # measurements_substrate = []
    # for data, init_substrate in zip(substrate, initial_substrates):
    #     measurement = {}
    #     measurement["initial_substrate"] = init_substrate
    #     reps = []
    #     for replicate in data:
    #         rep_dict = {}
    #         rep_dict["replicate"] = replicate.tolist()
    #         reps.append(rep_dict)
    #     measurement["data"] = reps

    #     measurements_substrate.append(measurement)

    # measurements_product = []
    # for data in product:
    #     measurement = {}
    #     measurement["initial_substrate"] = 0
    #     reps = []
    #     for replicate in data:
    #         rep_dict = {}
    #         rep_dict["replicate"] = replicate.tolist()
    #         reps.append(rep_dict)
    #     measurement["data"] = reps

    #     measurements_product.append(measurement)

    # data_dict = {
    #     "pH": pH,
    #     "name": f"ABTS oxidation pH {pH} and {temperature}°C",
    #     "date": str(
    #         datetime(
    #             *[int(x) for x in date.split("-")],
    #         )
    #     ),
    #     "time": [to_seconds(x) for x in time],
    #     "temperature": temperature,
    #     "s0": measurements_substrate,
    #     "s1": measurements_product,
    #     "s2": measurements_substrate_control,
    #     "s3": measurements_product_control,
    # }
    # return data_dict


def blockshaped(arr: np.ndarray, nrows: int, ncols: int):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size
    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    assert h % nrows == 0, f"{h} rows is not evenly divisible by {nrows}"
    assert w % ncols == 0, f"{w} cols is not evenly divisible by {ncols}"
    return (
        arr.reshape(h // nrows, nrows, -1, ncols)
        .swapaxes(1, 2)
        .reshape(-1, nrows, ncols)
    )


def read_photometer(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="delimiter", encoding="utf-16", engine="python")


def string_to_float(string: str) -> float:
    number = re.sub(r"[^0-9.]", "", string)
    if len(number) == 0:
        return float("nan")
    else:
        return float(number)


def to_seconds(string: str) -> List[int]:
    time = datetime.strptime(string, "%H:%M:%S").time()
    return time.hour * 3600 + time.minute * 60 + time.second


def _coordinates_to_id(x: int, y: int) -> str:
    return f"{chr(y + 65)}{x+1}"
