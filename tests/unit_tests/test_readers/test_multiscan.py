import numpy as np

from MTPHandler.model import Plate
from MTPHandler.readers import read_multiskan_spectrum_1500
from MTPHandler.units import C, minute

time = np.linspace(0, 15, 31).tolist()

ph = 6.9


def test_read_multiskan_spectrum_1500():
    # Arrange
    path = "docs/examples/data/multiskan_spectrum_1500.txt"
    print(time)

    # Act
    plate = read_multiskan_spectrum_1500(
        path=path,
        time=time,
        time_unit=minute,
        temperature=37.0,
        temperature_unit=C,
        ph=ph,
    )

    # Assert

    assert isinstance(plate, Plate)
    assert plate.temperatures == [37.0]
    assert plate.time_unit.name == minute.name
    assert plate.times == time
    assert len(plate.wells) == 60

    for well in plate.wells:
        if well.id == "A1":
            assert well.ph == ph
            assert well.x_pos == 0
            assert well.y_pos == 0
            measurment = well.measurements[0]
            assert measurment.wavelength == 340.0
            assert measurment.time_unit.name == minute.name
            assert measurment.time == time
            assert measurment.absorption[0] == 0.0440802853056135
            assert len(measurment.absorption) == len(time)
        if well.id == "A2":
            assert well.ph == ph
            assert well.x_pos == 1
            assert well.y_pos == 0
            measurment = well.measurements[0]
            assert measurment.wavelength == 340.0
            assert measurment.time_unit.name == minute.name
            assert measurment.time == time
            assert measurment.absorption[1] == 0.116455329304167
            assert len(measurment.absorption) == len(time)
