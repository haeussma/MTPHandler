import numpy as np
import pytest
from devtools import pprint

from MTPHandler import PlateManager
from MTPHandler.model import Plate
from MTPHandler.units import C, mmol, nm, s, ul

# create an artificaial plate with some wells and species


@pytest.fixture
def setup_handler() -> tuple[PlateManager, Plate]:
    plate = Plate(
        id="ID123",
        name="Test Plate",
        date_measured="2021-01-01 12:00:00",
        temperature_unit=C,
        temperatures=[37.0] * 12,
        time_unit=s,
        times=np.arange(0, 12, 1).tolist(),
    )
    w1 = plate.add_to_wells(
        id="A1",
        x_pos=0,
        y_pos=0,
        ph=7.0,
        volume=100,
        volume_unit=ul,
    )

    w2 = plate.add_to_wells(
        id="B1",
        x_pos=0,
        y_pos=1,
        ph=7.7,
        volume=200,
        volume_unit=ul,
    )

    w3 = plate.add_to_wells(
        id="A9",
        x_pos=8,
        y_pos=0,
        ph=8.0,
        volume=300,
        volume_unit=ul,
    )

    wells = [w1, w2, w3]

    for well in wells:
        well.add_to_measurements(
            wavelength=600,
            wavelength_unit=nm,
            absorption=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
            time_unit=s,
            time=np.arange(0, 12, 1).tolist(),
        )

    handler = PlateManager(plate=plate)
    return handler, plate


def test_add_well(setup_handler):
    handler, plate = setup_handler

    # Test adding a new well
    well = handler.plate.add_to_wells(
        id="C3", x_pos=2, y_pos=2, ph=7.5, volume=400, volume_unit=ul
    )
    assert well.id == "C3"
    assert well.x_pos == 2
    assert well.y_pos == 2
    assert well.ph == 7.5
    assert well.volume == 400
    assert well.volume_unit == ul
    assert len(plate.wells) == 4

    # test adding measurements to the well
    well.add_to_measurements(
        wavelength=600,
        wavelength_unit=nm,
        absorption=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
        time_unit=s,
        time=np.arange(0, 12, 1).tolist(),
    )
    assert len(well.measurements) == 1
    assert well.measurements[0].wavelength == 600
    assert well.measurements[0].wavelength_unit == nm
    assert well.measurements[0].time_unit == s
    assert well.measurements[0].time == np.arange(0, 12, 1).tolist()
    assert well.measurements[0].absorption == [
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.1,
        1.2,
    ]


def test_add_species(setup_handler):
    handler, plate = setup_handler

    # Test adding a new species
    species = handler.add_species(ld_id="LD123", id="ID123", name="Test Species")
    assert species.ld_id == "LD123"
    assert species.id == "ID123"
    assert species.name == "Test Species"
    assert len(plate.species) == 1

    # Test updating an existing species
    updated_species = handler.add_species(
        ld_id="LD123", id="ID123", name="Updated Species"
    )
    assert updated_species.name == "Updated Species"
    assert len(plate.species) == 1  # Still one species


def test_add_molecule(setup_handler):
    handler, plate = setup_handler

    # Test adding a new molecule
    molecule = handler.add_molecule(
        ld_id="LD456",
        id="ID456",
        name="Test Molecule",
        smiles="C",
        inchi_key="InChIKey",
    )
    assert molecule.ld_id == "LD456"
    assert molecule.id == "ID456"
    assert molecule.name == "Test Molecule"
    assert molecule.smiles == "C"
    assert molecule.inchi_key == "InChIKey"
    assert len(plate.species) == 1

    # Test updating an existing molecule
    updated_molecule = handler.add_molecule(
        ld_id="LD456",
        id="ID456",
        name="Updated Molecule",
        smiles="C",
        inchi_key="InChIKey",
    )
    assert updated_molecule.name == "Updated Molecule"
    assert len(plate.species) == 1  # Still one molecule


def test_add_protein(setup_handler):
    handler, plate = setup_handler

    # Test adding a new protein
    protein = handler.add_protein(
        ld_id="LD789",
        id="ID789",
        name="Test Protein",
        sequence="MVLTIYPDELVQIVSDKK",
        organism="E. coli",
        organism_tax_id=83333,
    )
    assert protein.ld_id == "LD789"
    assert protein.id == "ID789"
    assert protein.name == "Test Protein"
    assert protein.sequence == "MVLTIYPDELVQIVSDKK"
    assert protein.organism == "E. coli"
    assert protein.organism_tax_id == 83333
    assert len(plate.species) == 1

    # Test updating an existing protein
    updated_protein = handler.add_protein(
        ld_id="LD789",
        id="ID789",
        name="Updated Protein",
        sequence="MVLTIYPDELVQIVSDKK",
        organism="E. coli",
        organism_tax_id=83333,
    )
    assert updated_protein.name == "Updated Protein"
    assert len(plate.species) == 1  # Still one protein


# test assign logic for species


def test_assing_to_all_wells(setup_handler):
    handler, plate = setup_handler

    # Test assigning a species to all wells
    species = handler.add_species(ld_id="LD1", id="ID1", name="Test Species1")
    species2 = handler.add_species(ld_id="LD2", id="ID2", name="Test Species2")

    handler.assign_to_all(
        species=species, init_conc=1.2, conc_unit=mmol, contributes_to_signal=None
    )
    for well in plate.wells:
        assert well.measurements[0].blank_states[0].species_id == species.id
        assert well.measurements[0].blank_states[0].contributes_to_signal == True

        assert well.init_conditions[0].species_id == species.id
        assert well.init_conditions[0].init_conc == 1.2

    handler.assign_to_all(
        species=species2, init_conc=0, conc_unit=mmol, contributes_to_signal=None
    )

    for well in plate.wells:
        assert well.measurements[0].blank_states[1].species_id == species2.id
        assert well.measurements[0].blank_states[1].contributes_to_signal == False

        assert well.init_conditions[1].species_id == species2.id
        assert well.init_conditions[1].init_conc == 0


def test_assign_to_columns(setup_handler):
    handler, plate = setup_handler

    # Test assigning a species to a column
    species = handler.add_species(ld_id="LD1", id="ID1", name="Test Species1")
    species2 = handler.add_species(ld_id="LD2", id="ID2", name="Test Species2")

    handler.assign_to_columns(
        species=species,
        init_concs=[1, 2],
        conc_unit=mmol,
        contributes_to_signal=None,
        column_ids=[1],
    )

    for well in plate.wells:
        if well.x_pos == 0:
            assert well.measurements[0].blank_states[0].species_id == species.id
            assert well.measurements[0].blank_states[0].contributes_to_signal == True

            assert well.init_conditions[0].species_id == species.id
        elif well.id == "A1":
            assert well.init_conditions[0].init_conc == 1
        elif well.id == "A2":
            assert well.init_conditions[0].init_conc == 2

        else:
            assert len(well.measurements[0].blank_states) == 0
            assert len(well.init_conditions) == 0

    handler.assign_to_columns(
        species=species2,
        init_concs=[12],
        conc_unit=mmol,
        contributes_to_signal=None,
        column_ids=[9],
    )

    for well in plate.wells:
        if well.id == "A9":
            assert well.measurements[0].blank_states[0].species_id == species2.id
            assert well.measurements[0].blank_states[0].contributes_to_signal == True

            assert well.init_conditions[0].species_id == species2.id
            assert well.init_conditions[0].init_conc == 12

    # test number of wells in column and init concs does not match
    with pytest.raises(AssertionError):
        handler.assign_to_columns(
            species=species,
            init_concs=[1, 2],
            conc_unit=mmol,
            contributes_to_signal=None,
            column_ids=[2],
        )
