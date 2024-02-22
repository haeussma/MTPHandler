# MTPHandler - Tool for Microtiter Plate Data Handling

[![Tests](https://github.com/FAIRChemistry/MTPHandler/actions/workflows/tests.yml/badge.svg)](https://github.com/FAIRChemistry/MTPHandler/actions/workflows/tests.yml)

## 🛤 What is MTPHandler?

MTPHandler is a tool for managing and processing data from microtiter plates. Central to this tool is the `Plate` object, which provides comprehensive methods for manipulating chemical species within microtiter plates. This includes the addition and removal of species, assigning them to individual wells, and setting initial concentrations of reaction components.

Key Features of MTPHandler:

- __Parser functions__: Features a custom parser for different plate readers such as SpectraMax, Megellan, and MultiScan photometers, allowing for the mapping of raw TXT file data into the `Plate` data model. More photometers will be supported in the future.
- __Adaptive Data Processing__: Automatically adapts and blanks measurement data based on initial conditions set for each well.
- __Data Integration__: Incorporates additional data like pH and reaction temperature into the Plate object, which is not present in the photometer's TXT file.
- __Enhanced Object Definitions__: Utilizes Reactant and Protein objects, akin to those in EnzymeML Documents, complete with Systems Biology Ontology (SBO) annotations.
- __Versatile Experimental Applications__: Separately treats wells with and without protein species, facilitating both calibration (for creating standard curves) and reaction analysis.
- __EnzymeML Integration__: Maps well data to the EnzymeML data model using the specified conditions of the `Well` objects.

## ⚡️ Quick Start

Get started with CaliPytion by cloning this repository:

```Bash
git clone https://github.com/FAIRChemistry/MTPHandler/

```

Or install from PyPi:

```Bash
```

## 🔖 Example Code

coming soon

## ⚖️ License

Copyright (c) 2023 FAIR Chemistry

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
