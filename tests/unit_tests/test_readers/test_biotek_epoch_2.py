import pytest

from MTPHandler.model import Plate
from MTPHandler.readers import read_biotek
from MTPHandler.units import C, minute

ph = 6.9


def test_read_biotek_epoch_2():
    # Arrange
    path = "docs/examples/data/BioTek_Epoch2.xlsx"
    # Act
    plate = read_biotek(
        path=path,
        ph=ph,
    )

    # Assert
    assert isinstance(plate, Plate)
    assert plate.temperature_unit.name == C.name
    assert len(plate.wells) == 36

    for well in plate.wells:
        if well.id == "B2":
            assert len(well.measurements) == 1
            assert well.ph == ph
            assert well.x_pos == 1
            assert well.y_pos == 1
            measurment = well.measurements[0]
            assert measurment.wavelength == 630
            assert measurment.time_unit.name == minute.name
            assert measurment.absorption[4] == pytest.approx(0.1, rel=1e-2)
            assert len(measurment.absorption) == 353

        if well.id == "D9":
            assert well.measurements[0].absorption[-1] == pytest.approx(0.226, rel=1e-2)
