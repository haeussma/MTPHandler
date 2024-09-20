import pytest

from MTPHandler.model import Plate
from MTPHandler.readers import read_spectra_max_190
from MTPHandler.units import C, minute

ph = 6.9


def test_spectramax_190():
    # Arrange
    path1 = "docs/examples/data/Spectramax190 molecular Devices.txt"
    path2 = "docs/examples/data/spectra_max_190.txt"

    # Act
    plate1 = read_spectra_max_190(
        path=path1,
        ph=ph,
    )

    plate2 = read_spectra_max_190(
        path=path2,
        ph=ph,
    )

    # Assert plate1
    assert isinstance(plate1, Plate)
    assert plate1.temperature_unit.name == C.name
    assert plate1.temperatures[0] == 21.6
    assert len(plate1.wells) == 96

    for well in plate1.wells:
        if well.id == "A3":
            assert len(well.measurements) == 1
            assert well.ph == ph
            assert well.x_pos == 2
            assert well.y_pos == 0
            measurment = well.measurements[0]
            assert measurment.wavelength == 500.0
            assert measurment.time_unit.name == minute.name
            assert measurment.absorption[4] == pytest.approx(0.099, rel=1e-2)
            assert len(measurment.absorption) == 9

    # assert plate2
    assert isinstance(plate2, Plate)
    assert plate2.temperature_unit.name == C.name
    assert plate2.temperatures[0] == 37.0
    assert len(plate2.wells) == 96

    for well in plate2.wells:
        if well.id == "A3":
            assert len(well.measurements) == 1
            assert well.ph == ph
            assert well.x_pos == 2
            assert well.y_pos == 0
            measurment = well.measurements[0]
            assert measurment.wavelength == 412.0
            assert measurment.time_unit.name == minute.name
            assert measurment.absorption[4] == pytest.approx(0.039, rel=1e-2)
            assert len(measurment.absorption) == 721
