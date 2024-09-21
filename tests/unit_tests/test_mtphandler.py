import numpy as np
import pytest

<<<<<<< HEAD
from mtphandler import PlateManager
from mtphandler.model import Plate
from mtphandler.molecule import Molecule
from mtphandler.units import C, mmol, s, ul


# create an artificaial plate with some wells and species
@pytest.fixture
def get_molecule():
    return Molecule(
        id="m1",
        pubchem_cid=12345,
        name="Test Molecule",
        constant=False,
    )


@pytest.fixture
def get_other_molecule():
    return Molecule(
        id="m2",
        pubchem_cid=123456,
        name="Test Molecule 2",
        constant=False,
    )
=======
from MTPHandler import PlateManager
from MTPHandler.model import Plate
from MTPHandler.units import C, mmol, nm, s, ul

# create an artificaial plate with some wells and species
>>>>>>> main


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
<<<<<<< HEAD
=======
            wavelength_unit=nm,
>>>>>>> main
            absorption=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
            time_unit=s,
            time=np.arange(0, 12, 1).tolist(),
        )

<<<<<<< HEAD
    handler = PlateManager(plate=plate, name="Test Plate")
=======
    handler = PlateManager(plate=plate)
>>>>>>> main
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
<<<<<<< HEAD
=======
        wavelength_unit=nm,
>>>>>>> main
        absorption=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
        time_unit=s,
        time=np.arange(0, 12, 1).tolist(),
    )
    assert len(well.measurements) == 1
    assert well.measurements[0].wavelength == 600
<<<<<<< HEAD
=======
    assert well.measurements[0].wavelength_unit == nm
>>>>>>> main
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


<<<<<<< HEAD
def test_assing_to_all_wells(setup_handler, get_molecule, get_other_molecule):
    handler, plate = setup_handler

    mol1 = get_molecule
    mol2 = get_other_molecule

    # Test assigning a species to all wells
    handler.add_molecule(molecule=mol1)
    handler.add_molecule(molecule=mol2)

    handler._assign_to_all(
        species=mol1,
        init_conc=1.2,
        conc_unit=mmol,
        contributes_to_signal=None,
        silent=False,
    )
    for well in plate.wells:
        assert well.measurements[0].blank_states[0].species_id == mol1.id
        assert well.measurements[0].blank_states[0].contributes_to_signal

        assert well.init_conditions[0].species_id == mol1.id
        assert well.init_conditions[0].init_conc == 1.2

    handler._assign_to_all(
        species=mol2,
        init_conc=0,
        conc_unit=mmol,
        contributes_to_signal=None,
        silent=False,
    )

    for well in plate.wells:
        assert well.measurements[0].blank_states[1].species_id == mol2.id
        assert not well.measurements[0].blank_states[1].contributes_to_signal

        assert well.init_conditions[1].species_id == mol2.id
        assert well.init_conditions[1].init_conc == 0


def test_assign_to_columns(setup_handler, get_molecule, get_other_molecule):
    handler, plate = setup_handler

    mol1 = get_molecule
    mol2 = get_other_molecule

    # Test assigning a species to all wells
    handler.add_molecule(molecule=mol1)
    handler.add_molecule(molecule=mol2)

    handler._assign_to_columns(
        species=mol1,
=======
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
>>>>>>> main
        init_concs=[1, 2],
        conc_unit=mmol,
        contributes_to_signal=None,
        column_ids=[1],
<<<<<<< HEAD
        silent=False,
=======
>>>>>>> main
    )

    for well in plate.wells:
        if well.x_pos == 0:
<<<<<<< HEAD
            assert well.measurements[0].blank_states[0].species_id == mol1.id
            assert well.measurements[0].blank_states[0].contributes_to_signal

            assert well.init_conditions[0].species_id == mol1.id
=======
            assert well.measurements[0].blank_states[0].species_id == species.id
            assert well.measurements[0].blank_states[0].contributes_to_signal == True

            assert well.init_conditions[0].species_id == species.id
>>>>>>> main
        elif well.id == "A1":
            assert well.init_conditions[0].init_conc == 1
        elif well.id == "A2":
            assert well.init_conditions[0].init_conc == 2

        else:
            assert len(well.measurements[0].blank_states) == 0
            assert len(well.init_conditions) == 0

<<<<<<< HEAD
    handler._assign_to_columns(
        species=mol2,
=======
    handler.assign_to_columns(
        species=species2,
>>>>>>> main
        init_concs=[12],
        conc_unit=mmol,
        contributes_to_signal=None,
        column_ids=[9],
<<<<<<< HEAD
        silent=False,
=======
>>>>>>> main
    )

    for well in plate.wells:
        if well.id == "A9":
<<<<<<< HEAD
            assert well.measurements[0].blank_states[0].species_id == mol2.id
            assert well.measurements[0].blank_states[0].contributes_to_signal

            assert well.init_conditions[0].species_id == mol2.id
=======
            assert well.measurements[0].blank_states[0].species_id == species2.id
            assert well.measurements[0].blank_states[0].contributes_to_signal == True

            assert well.init_conditions[0].species_id == species2.id
>>>>>>> main
            assert well.init_conditions[0].init_conc == 12

    # test number of wells in column and init concs does not match
    with pytest.raises(AssertionError):
<<<<<<< HEAD
        handler._assign_to_columns(
            species=mol1,
=======
        handler.assign_to_columns(
            species=species,
>>>>>>> main
            init_concs=[1, 2],
            conc_unit=mmol,
            contributes_to_signal=None,
            column_ids=[2],
<<<<<<< HEAD
            silent=False,
        )


def test_add_molecule_without_constant(setup_handler):
    plate_manager, plate = setup_handler

    molecule = Molecule(
        id="m1", pubchem_cid=12345, name="Test Molecule", constant=False
    )

    plate_manager.add_molecule(molecule)

    assert len(plate_manager.molecules) == 1
    assert plate_manager.molecules[0].id == "m1"
    assert plate_manager.molecules[0].constant == False


def test_update_existing_molecule(setup_handler):
    plate_manager, plate = setup_handler

    molecule1 = Molecule(
        id="m1", pubchem_cid=12345, name="Test Molecule", constant=False
    )
    molecule2 = Molecule(
        id="m1", pubchem_cid=12345, name="Updated Molecule", constant=True
    )

    plate_manager.add_molecule(molecule1)
    plate_manager.add_molecule(molecule2)

    assert len(plate_manager.molecules) == 1
    assert plate_manager.molecules[0].id == "m1"
    assert plate_manager.molecules[0].name == "Updated Molecule"
    assert plate_manager.molecules[0].constant == True


def test_define_molecule_with_valid_data(setup_handler):
    plate_manager, plate = setup_handler

    molecule = plate_manager.define_molecule(
        id="m1", pubchem_cid=12345, name="Test Molecule", constant=False
    )

    assert molecule.id == "m1"
    assert molecule.pubchem_cid == 12345
    assert molecule.name == "Test Molecule"
    assert molecule.constant is False


def test_define_molecule_with_non_integer_pubchem_cid(setup_handler):
    plate_manager, plate = setup_handler

    with pytest.raises(ValueError, match="PubChem CID must be an integer."):
        plate_manager.define_molecule(
            id="m1",
            pubchem_cid="not_an_integer",
            name="Test Molecule",
            constant=False,
        )


def test_define_molecule_without_name(setup_handler):
    # Initialize the class object
    plate_manager, plate = setup_handler

    # Define a molecule with valid id and pubchem_cid, but without name
    molecule_id = "m1"
    pubchem_cid = 12345

    mol = plate_manager.define_molecule(id=molecule_id, pubchem_cid=pubchem_cid)

    assert len(mol.name) > 0


def test_append_molecule_if_not_exist(setup_handler):
    # Initialize the PlateManager object
    plate_manager, plate = setup_handler

    # Define a molecule with PubChem CID not available
    molecule = plate_manager.define_molecule(
        id="m1", pubchem_cid=-1, name="Test Molecule", constant=False
    )

    # Check if the molecule is correctly appended to the list
    assert len(plate_manager.molecules) == 1
    assert plate_manager.molecules[0].id == "m1"
    assert plate_manager.molecules[0].name == "Test Molecule"
    assert plate_manager.molecules[0].constant == False


def test_defining_already_defined_molecule_updates_molecule(setup_handler):
    # Initialize the PlateManager object
    plate_manager, plate = setup_handler

    # Define a molecule with PubChem CID not available
    molecule = plate_manager.define_molecule(
        id="m1", pubchem_cid=-1, name="Test Molecule", constant=False
    )

    # Check if the molecule is correctly appended to the list
    assert len(plate_manager.molecules) == 1
    assert plate_manager.molecules[0].id == "m1"
    assert plate_manager.molecules[0].name == "Test Molecule"
    assert plate_manager.molecules[0].constant == False

    # Define the molecule again with different properties
    molecule = plate_manager.define_molecule(
        id="m1", pubchem_cid=12345, name="New Molecule", constant=True
    )

    # Check if the molecule is updated
    assert len(plate_manager.molecules) == 1
    assert plate_manager.molecules[0].id == "m1"
    assert plate_manager.molecules[0].name == "New Molecule"
    assert plate_manager.molecules[0].constant == True


def test_define_molecule_pubchem_cid_minus_one_no_name(setup_handler):
    # Initialize the class object
    plate_manager, plate = setup_handler

    # Define molecule with pubchem_cid=-1 and no name provided
    with pytest.raises(ValueError):
        molecule = plate_manager.define_molecule(id="test_id", pubchem_cid=-1)


def test_molecule_string_pubchem_cid(setup_handler):
    # Initialize the class object
    plate_manager, plate = setup_handler

    # Define molecule with invalid pubchem_cid
    with pytest.raises(ValueError, match="PubChem CID must be an integer."):
        molecule = plate_manager.define_molecule(
            id="test_id", pubchem_cid="invalid_pubchem_cid"
        )


def test_molecule_cid_not_exist(setup_handler):
    # Initialize the class object
    plate_manager, plate = setup_handler

    # Define molecule with pubchem_cid=-1 and no name provided
    with pytest.raises(
        ValueError, match="Failed to retrieve molecule name from PubChem"
    ):
        molecule = plate_manager.define_molecule(
            id="test_id", pubchem_cid=12312323235235235235
=======
>>>>>>> main
        )
