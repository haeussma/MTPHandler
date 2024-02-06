# MTPHandler

## ⚡️ Quick Start

To install MTPHandler, run the following command:

```bash
pip install git+https://github.com/pi/pa/po.git
```

```python
from MTPHandler import Plate
from MTPHandler.parsers import SpectraMaxParser

# Create a Plate object
plate = Plate.parse(SpectraMaxParser, "path/to/your/data.txt", ph=7.3)

# Visualize the plate
plate.visualize()
```
