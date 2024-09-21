<<<<<<< HEAD
from mtphandler.model import Plate
from mtphandler.readers import read_tekan_magellan
from mtphandler.units import second
=======
from MTPHandler.model import Plate
from MTPHandler.readers import read_tekan_magellan
from MTPHandler.units import nm, s
>>>>>>> main


def test_magellan():
    # Arrange
<<<<<<< HEAD
    path = "docs/examples/data/magellan.xlsx"
=======
    path = "tests/data/magellan.xlsx"
>>>>>>> main
    ph = 7.0
    wavelength = 450.0

    # Act
    plate = read_tekan_magellan(
        path=path,
        wavelength=wavelength,
        ph=ph,
    )

    # Assert
<<<<<<< HEAD

    assert isinstance(plate, Plate)
    assert plate.temperatures == [37.1, 37.2, 37.1, 37.2, 36.9, 37.3, 37.2]
    assert plate.time_unit.name == second.name
=======
    assert isinstance(plate, Plate)
    assert plate.temperatures == [37.1, 37.2, 37.1, 37.2, 36.9, 37.3, 37.2]
    assert plate.time_unit == s
>>>>>>> main
    assert plate.times == [0, 900.0, 1802.0, 2700.0, 3600.0, 4502.0, 5400.0]
    assert len(plate.wells) == 62

    for well in plate.wells:
        if well.id == "E8":
            assert well.ph == ph
            assert well.x_pos == 7
            assert well.y_pos == 4
            measurment = well.measurements[0]
            assert measurment.wavelength == 450.0
<<<<<<< HEAD
            assert measurment.time_unit.name == second.name
=======
            assert measurment.wavelength_unit == nm
            assert measurment.time_unit == s
>>>>>>> main
            assert measurment.time == [0, 900.0, 1802.0, 2700.0, 3600.0, 4502.0, 5400.0]
            assert measurment.absorption[0] == 2.6871
            assert len(measurment.absorption) == 7
