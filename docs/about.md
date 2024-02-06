# What is MTPHandler

MTPHandler is a tool for managing and processing data from microtiter plates. Central to this tool is the `Plate` object, which provides comprehensive methods for manipulating chemical species within microtiter plates. This includes the addition and removal of species, assigning them to individual wells, and setting initial concentrations of reaction components.

Key Features of MTPHandler:

- __Parser functions__: Features a custom parser for different plate readers such as SpectraMax, Megellan, and MultiScan photometers, allowing for the mapping of raw TXT file data into the `Plate` data model. More photometers will be supported in the future.
- __Adaptive Data Processing__: Automatically adapts and blanks measurement data based on initial conditions set for each well.
- __Data Integration__: Incorporates additional data like pH and reaction temperature into the Plate object, which is not present in the photometer's TXT file.
- __Enhanced Object Definitions__: Utilizes Reactant and Protein objects, akin to those in EnzymeML Documents, complete with Systems Biology Ontology (SBO) annotations.
- __Versatile Experimental Applications__: Separately treats wells with and without protein species, facilitating both calibration (for creating standard curves) and reaction analysis.
- __EnzymeML Integration__: Maps well data to the EnzymeML data model using the specified conditions of the `Well` objects.

## 🏁 Currently implemented photometers:

- [x]Spectramax
- [x] Megellan
- [x] MultiScan