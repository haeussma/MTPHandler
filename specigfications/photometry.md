# Spectrophotometry data model

## Modules

### Plate

- n_rows
    - Type: integer
    - Description: Number of rows on the plate
- n_columns
    - Type: integer
    - Description: Number of columns on the plate
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
- init_conc
    - Type: float
    - Description: Initial concentration of the species
- conc_unit
    - Type: str
    - Description: Concentration unit
- x_position
    - Type: integer
    - Description: X position of the well on the plate
- y_position
    - Type: integer
    - Description: Y position of the well on the plate
- species_id
    - Type: str
    - Description: ID of the species
- wavelegth
    - Type: int
    - Description: Wavelength of the measurement
