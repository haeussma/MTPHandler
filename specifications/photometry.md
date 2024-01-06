# Data model for microtiter plate photometry

## Modules

### Plate

- __n_rows__
    - Type: integer
    - Description: Number of rows on the plate
- __n_columns__
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
- __temperatures__
    - Type: float
    - Description: Thermostat temperature
    - Multiple: True
- __temperature_unit__
    - Type: str
    - Description: Unit of the temperature
- max_volume
    - Type: float
    - Description: Maximum volume of the wells
- max_volume_unit
    - Type: str
    - Description: Unit of the maximum volume
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
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@AbstractSpecies
    - Description: List of species present in wells of the plate
    - Multiple: True

### Well

- __ph__
    - Type: float
    - Description: pH of the reaction
- __x_position__
    - Type: integer
    - Description: X position of the well on the plate
- __y_position__
    - Type: integer
    - Description: Y position of the well on the plate
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

### PhotometricMeasurement

- __wavelength__
    - Type: float
    - Description: Wavelength of the measurement
- __wavelength_unit__
    - Type: str
    - Description: Unit of the wavelength
- __absorptions__
    - Type: float
    - Description: Absorption of the species
    - Multiple: True
- blank_states
    - Type: BlankState
    - Description: List of blank states, referring to the blank state of the species of the well
    - Multiple: True

### InitCondition

- __species_id__
    - Type: @AbstractSpecies.id
    - Description: Reference to species
- __init_conc__
    - Type: float
    - Description: Initial concentration of the species
- __conc_unit__
    - Type: str
    - Description: Concentration unit

### BlankState

- __species_id__
    - Type: @AbstractSpecies.id
    - Description: Reference to species
- __contributes_to_signal__
    - Type: bool
    - Description: Whether the species' absorption contributes to the absorption signal
    - Default: True
