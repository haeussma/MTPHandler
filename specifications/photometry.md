# Data model for microtiter plate photometry

## Modules

### Plate

Description of a microtiter plate consisting of wells.

- __n_rows__
    - Type: integer
    - Description: Number of rows on the plate
- __n_cols__
    - Type: integer
    - Description: Number of columns on the plate
- date_measured
    - Type: datetime
    - Description: Date and time when the plate was measured
- times
    - Type: float[]
    - Description: Time points of the measurement, corresponding to temperature measurements
- time_unit
    - Type: Unit
    - Description: Unit of the time
- temperatures
    - Type: float[]
    - Description: Thermostat temperature
- temperature_unit
    - Type: Unit
    - Description: Unit of the temperature
- wells
    - Type: Well[]
    - Description: List of wells on the plate
- species
    - Type: Reactant[], Protein[]
    - Description: List of species present in wells of the plate

### Protein

Description of a protein species that might be present in the wells of the plate.

- name
    - Type: str
    - Description: Name of the species
- sequence
    - Type: str
    - Description: Amino acid sequence of the protein
- organism
    - Type: str
    - Description: Organism the protein originates from
- organism_tax_id
    - Type: Identifier
    - Description: NCBI taxonomy ID of the organism
- references
    - Type: Identifier[]
    - Description: List of references to the protein

### Reactant

Description of a chemical species that might be present in the wells of the plate.

- name
    - Type: str
    - Description: Name of the species
- smiles
    - Type: str
    - Description: SMILES representation of the species
- inchi
    - Type: str
    - Description: InChI representation of the species
- references
    - Type: Identifier[]
    - Description: List of references to the Reactant

### Well

Description of a well on the plate.

- __x_pos__
    - Type: integer
    - Description: X position of the well on the plate
- __y_pos__
    - Type: integer
    - Description: Y position of the well on the plate
- ph
    - Type: float
    - Description: pH of the reaction
- init_conditions
    - Type: InitCondition[]
    - Description: List of initial conditions of different species
- measurements
    - Type: PhotometricMeasurement[]
    - Description: List of photometric measurements
- volume
    - Type: float
    - Description: Volume of the reaction
- volume_unit
    - Type: Unit
    - Description: Unit of the volume

### PhotometricMeasurement

Description of a photometric measurement of a single well on the plate.

- __wavelength__
    - Type: float
    - Description: Wavelength of the measurement
- __wavelength_unit__
    - Type: Unit
    - Description: Unit of the wavelength
- absorption
    - Type: float[]
    - Description: Absorption of the species
- time
    - Type: float
    - Description: Time of the measurement
- time_unit
    - Type: Unit
    - Description: Unit of the time
- blank_states
    - Type: BlankState[]
    - Description: List of blank states, referring to the blank state of the species of the well

### InitCondition

Description of the initial condition of a species in a well.

- __species_id__
    - Type: str
    - Description: Reference to species
- __init_conc__
    - Type: float
    - Description: Initial concentration of the species
- __conc_unit__
    - Type: Unit
    - Description: Concentration unit

### BlankState

Describes if the respective species contributes to the absorption signal.

- __species_id__
    - Type: str
    - Description: Reference to species
- __contributes_to_signal__
    - Type: bool
    - Description: Whether the species' absorption contributes to the absorption signal
    - Default: True
