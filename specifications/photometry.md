# Datamodel to handel microtiter plate photometry data

## Modules

### Plate

- n_rows
    - Type: integer
    - Description: Number of rows on the plate
- n_columns
    - Type: integer
    - Description: Number of columns on the plate
- date_measured
    - Type: datetime
    - Description: Date and time when the plate was measured
- times
    - Type: float
    - Description: Time points of the measurement, corresponding to temperature measurements
    - Multiple: True
- time_unit
    - Type: str
    - Description: Unit of the time
- temperatures
    - Type: float
    - Description: Thermostat temperature
    - Multiple: True
- temperature_unit
    - Type: str
    - Description: Unit of the temperature
- max_volume
    - Type: float
    - Description: Maximum volume of the wells
- max_volume_unit
    - Type: str
    - Description: Unit of the maximum volume
- ph
    - Type: float
    - Description: pH of the reaction
- wells
    - Type: Well
    - Description: List of wells on the plate
    - Multiple: True
- measured_wavelengths
    - Type: float
    - Description: Measured wavelengths
    - Multiple: True
- wavelength_unit
    - Type: str
    - Description: Unit of the wavelength
- species
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@AbstractSpecies, https://github.com/EnzymeML/enzymeml-specifications.git@Protein, https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: List of species present in wells of the plate

### Well

- init_conditions
    - Type: InitCondition
    - Multiple: True
    - Description: List of initial conditions of different species
- measurements
    - Type: PhotometricMeasurement
    - Multiple: True
    - Description: List of photometric measurements
- volume
    - Type: float
    - Description: Volume of the reaction
- volume_unit
    - Type: string
    - Description: Unit of the volume
- x_position
    - Type: integer
    - Description: X position of the well on the plate
- y_position
    - Type: integer
    - Description: Y position of the well on the plate
- wavelength
    - Type: int
    - Description: Wavelength of the measurement

### PhotometricMeasurement

- wavelength
    - Type: float
    - Description: Wavelength of the measurement
- wavelength_unit
    - Type: str
    - Description: Unit of the wavelength
- time
    - Type: float
    - Description: Time of the measurement
- time_unit
    - Type: str
    - Description: Unit of the time
- absorptions
    - Type: float
    - Description: Absorption of the species
    - Multiple: True

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