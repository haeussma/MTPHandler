---
repo: https://github.com/FAIRChemistry/MTPHandler
---

# Data model for microtiter plate photometry

## Modules

### Plate

Description of a microtiter plate consisting of wells.

- id
    - Type: string
    - Description: Identifier of the plate
- name
    - Type: string
    - Description: Arbitrary name of the plate.
- wells
    - Type: Well[]
    - Description: List of wells on the plate.
- date_measured
    - Type: string
    - Description: Date and time when the plate was measured.
- __temperatures__
    - Type: float[]
    - Description: Thermostat temperature
- __temperature_unit__
    - Type: UnitDefinition
    - Description: Unit of the temperature
- times
    - Type: float[]
    - Description: Time points of the measurement, corresponding to temperature measurements.
- time_unit
    - Type: UnitDefinition
    - Description: Unit of the time

### Well

Description of a well on the plate.

- __id__
    - Type: string
    - Description: Identifier of the well
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
    - Type: UnitDefinition
    - Description: Unit of the volume

### PhotometricMeasurement

Description of a photometric measurement of a single well and wavelength on the plate.

- __wavelength__
    - Type: float
    - Description: Wavelength of the measurement
- absorption
    - Type: float[]
    - Description: Absorption of the species
- time
    - Type: float[]
    - Description: Time of the measurement
- time_unit
    - Type: UnitDefinition
    - Description: Unit of the time
- blank_states
    - Type: BlankState[]
    - Description: List of blank states, referring to the blank state of the species of the well

### InitCondition

Description of the initial condition of a species in a well.

- __species_id__
    - Type: string
    - Description: Reference to species
- __init_conc__
    - Type: float
    - Description: Initial concentration of the species
- __conc_unit__
    - Type: UnitDefinition
    - Description: Concentration unit

### BlankState

Describes if the respective species contributes to the absorption signal.

- __species_id__
    - Type: string
    - Description: Reference to species
- __contributes_to_signal__
    - Type: boolean
    - Description: Whether the species' absorption contributes to the absorption signal
    - Default: True
