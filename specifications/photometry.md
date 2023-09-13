# Datamodel to handel microtiter plate photometry data

## Modules

### Plate

- n_rows
    - Type: integer
    - Description: Number of rows on the plate
- n_columns
    - Type: integer
    - Description: Number of columns on the plate
- created
    - Type: datetime
    - Description: Date and time when the plate was measured
- temperature
    - Type: float
    - Description: Thermostat temperature
- temperature_unit
    - Type: str
    - Description: Unit of the temperature
- ph
    - Type: float
    - Description: pH of the reaction
- wells
    - Type: Well
    - Description: List of wells on the plate
    - Multiple: True
- measured_wavelengths
    - Type: int
    - Description: Measured wavelengths in nm
    - Multiple: True
- species
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@AbstractSpecies, https://github.com/EnzymeML/enzymeml-specifications.git@Protein, https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: List of species present in wells of the plate

### Well

- absorption
    - Type: float
    - Multiple: True
    - Description: Absorption of the species
- time
    - Type: float
    - Multiple: True
    - Description: Time of the measurement
- time_unit
    - Type: str
    - Description: Unit of the time
- reaction_volume
    - Type: float
    - Description: Volume of the reaction
- volume_unit
    - Type: string
    - Description: Unit of the volume
- init_conditions
    - Type: InitCondition
    - Multiple: True
    - Description: List of initial conditions of different species
- x_position
    - Type: integer
    - Description: X position of the well on the plate
- y_position
    - Type: integer
    - Description: Y position of the well on the plate
- wavelength
    - Type: int
    - Description: Wavelength of the measurement

### InitCondition

- species_id
    - Type: @AbstractSpecies.id
    - Description: Reference to species
- init_conc
    - Type: float
    - Description: Initial concentration of the species
- conc_unit
    - Type: str
    - Description: Concentration unit
- __was_blanked__
    - Type: bool
    - Description: Whether the species' absorption contribution was subtracted from the absorption signal
    - Default: False

## Enumerations

### SpeciesType

```python
BUFFER = "buffer"
ENZYME = "enzyme"
SUBSTRATE = "substrate"
PRODUCT = "product"
INHIBITOR = "inhibitor"
```

