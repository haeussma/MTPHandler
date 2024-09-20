---
hide:
    - navigation
---

# `Plate`

!!! quote "Graph"
    ```mermaid
    flowchart TB
        plate(Plate)
        well(Well)
        photometricmeasurement(PhotometricMeasurement)
        initcondition(InitCondition)
        blankstate(BlankState)
        unitdefinition(UnitDefinition)
        baseunit(BaseUnit)
        unittype(UnitType)
        plate(Plate) --> well(Well)
        plate(Plate) --> unitdefinition(UnitDefinition)
        plate(Plate) --> unitdefinition(UnitDefinition)
        well(Well) --> initcondition(InitCondition)
        well(Well) --> photometricmeasurement(PhotometricMeasurement)
        well(Well) --> unitdefinition(UnitDefinition)
        photometricmeasurement(PhotometricMeasurement) --> unitdefinition(UnitDefinition)
        photometricmeasurement(PhotometricMeasurement) --> blankstate(BlankState)
        initcondition(InitCondition) --> unitdefinition(UnitDefinition)
        unitdefinition(UnitDefinition) --> baseunit(BaseUnit)
        baseunit(BaseUnit) --> unittype(UnitType)

        click plate "#plate" "Go to Plate"
        click well "#well" "Go to Well"
        click photometricmeasurement "#photometricmeasurement" "Go to PhotometricMeasurement"
        click initcondition "#initcondition" "Go to InitCondition"
        click blankstate "#blankstate" "Go to BlankState"
        click unitdefinition "#unitdefinition" "Go to UnitDefinition"
        click baseunit "#baseunit" "Go to BaseUnit"
        click unittype "#unittype" "Go to UnitType"
    ```


## Types


### Plate
Description of a microtiter plate consisting of wells.

__id__ `string`

- Identifier of the plate


__name__ `string`

- Arbitrary name of the plate.


__wells__ [`list[Well]`](#well)

- List of wells on the plate.


__date_measured__ `string`

- Date and time when the plate was measured.


__temperatures__* `list[float]`

- Thermostat temperature


__temperature_unit__* [`UnitDefinition`](#unitdefinition)

- Unit of the temperature


__times__ `list[float]`

- Time points of the measurement, corresponding to temperature measurements.


__time_unit__ [`UnitDefinition`](#unitdefinition)

- Unit of the time


------

### Well
Description of a well on the plate.

__id__* `string`

- Identifier of the well


__x_pos__* `integer`

- X position of the well on the plate


__y_pos__* `integer`

- Y position of the well on the plate


__ph__ `float`

- pH of the reaction


__init_conditions__ [`list[InitCondition]`](#initcondition)

- List of initial conditions of different species


__measurements__ [`list[PhotometricMeasurement]`](#photometricmeasurement)

- List of photometric measurements


__volume__ `float`

- Volume of the reaction


__volume_unit__ [`UnitDefinition`](#unitdefinition)

- Unit of the volume


------

### PhotometricMeasurement
Description of a photometric measurement of a single well and wavelength on the plate.

__wavelength__* `float`

- Wavelength of the measurement


__absorption__* `list[float]`

- Absorption of the species


__time__* `list[float]`

- Time of the measurement


__time_unit__* [`UnitDefinition`](#unitdefinition)

- Unit of the time


__blank_states__* [`list[BlankState]`](#blankstate)

- List of blank states, referring to the blank state of the species of the well


------

### InitCondition
Description of the initial condition of a species in a well.

__species_id__* `string`

- Reference to species


__init_conc__* `float`

- Initial concentration of the species


__conc_unit__* [`UnitDefinition`](#unitdefinition)

- Concentration unit


------

### BlankState
Describes if the respective species contributes to the absorption signal.

__species_id__* `string`

- Reference to species


__contributes_to_signal__* `boolean`

- Whether the species' absorption contributes to the absorption signal

- `Default`: true

------

### UnitDefinition
Represents a unit definition that is based on the SI unit system.

__id__ `string`

- Unique identifier of the unit definition.


__name__ `string`

- Common name of the unit definition.


__base_units__ [`list[BaseUnit]`](#baseunit)

- Base units that define the unit.


------

### BaseUnit
Represents a base unit in the unit definition.

__kind__* [`UnitType`](#unittype)

- Kind of the base unit (e.g., meter, kilogram, second).


__exponent__* `integer`

- Exponent of the base unit in the unit definition.


__multiplier__ `float`

- Multiplier of the base unit in the unit definition.


__scale__ `float`

- Scale of the base unit in the unit definition.


## Enumerations

### UnitType

| Alias | Value |
|-------|-------|
| `AMPERE` | ampere |
| `AVOGADRO` | avogadro |
| `BECQUEREL` | becquerel |
| `CANDELA` | candela |
| `CELSIUS` | celsius |
| `COULOMB` | coulomb |
| `DIMENSIONLESS` | dimensionless |
| `FARAD` | farad |
| `GRAM` | gram |
| `GRAY` | gray |
| `HENRY` | henry |
| `HERTZ` | hertz |
| `ITEM` | item |
| `JOULE` | joule |
| `KATAL` | katal |
| `KELVIN` | kelvin |
| `KILOGRAM` | kilogram |
| `LITRE` | litre |
| `LUMEN` | lumen |
| `LUX` | lux |
| `METRE` | metre |
| `MOLE` | mole |
| `NEWTON` | newton |
| `OHM` | ohm |
| `PASCAL` | pascal |
| `RADIAN` | radian |
| `SECOND` | second |
| `SIEMENS` | siemens |
| `SIEVERT` | sievert |
| `STERADIAN` | steradian |
| `TESLA` | tesla |
| `VOLT` | volt |
| `WATT` | watt |
| `WEBER` | weber |