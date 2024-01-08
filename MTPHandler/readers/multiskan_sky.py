import re
import pandas as pd


def read_multiskan_sky(
    cls: "Plate",
    path: str,
    ph: float,
    time_unit: str = "s",
):
    # Extract meta data of measurement
    df = pd.read_excel(path, header=None, sheet_name="Absorbance results")
    created = df.iloc[1, 0]

    temperature_cell = pd.read_excel(path, sheet_name="General information").iloc[6, 3]
    temperature = float(re.findall(r"\d+", temperature_cell)[0])

    # Check if there are multiple sheets with data
    data_sheets = []
    sheets = pd.ExcelFile(path).sheet_names
    for sheet in sheets:
        if "Raw data" in sheet:
            data_sheets.append(sheet)

    # Extract wavelengths
    wavelengths = []
    for sheet_name in data_sheets:
        wavelengths.append(int(re.findall(r"\d+", sheet_name)[0]))

    # Create plate
    plate = cls(
        measured_wavelengths=wavelengths,
        temperature=temperature,
        date_measured=created,
        temperature_unit="C",
        time_unit=time_unit,
        n_rows=8,
        n_columns=12,
    )

    # Extract data and add it to 'Well' objects
    for sheet_name in data_sheets:
        wavelength = int(re.findall(r"\d+", sheet_name)[0])
        df = pd.read_excel(path, header=1, sheet_name=sheet_name, skiprows=3)

        for _, row in df.iterrows():
            well_id = row.values[0].replace(" ", "")
            wavelength = row.values[1]
            absorption = row.values[2]
            time = row.values[3]

            # Add well to plate if it does not exist
            try:
                well = plate.get_well(well_id)
            except ValueError:
                x_pos, y_pos = id_to_xy(well_id)
                well = plate.add_to_wells(
                    id=well_id,
                    ph=ph,
                    x_position=x_pos,
                    y_position=y_pos,
                )

            # Add measurement to well
            try:
                measurement = well.get_measurement(wavelength)
            except ValueError:
                measurement = well.add_to_measurements(
                    wavelength=wavelength,
                    wavelength_unit="nm",
                )

            measurement.absorptions.append(absorption)
            measurement.times.append(time)

    return plate


def id_to_xy(well_id: str):
    return int(well_id[1:]) - 1, ord(well_id[0].upper()) - 65
