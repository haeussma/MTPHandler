import numpy as np

from MTPHandler.dataclasses import Plate
from MTPHandler.units import celsius, minute, mM

times = np.linspace(0, 10, 11)


p = Plate(
    id="MTP_001",
    name="Enzyme Kinetics",
    temperatures=[25] * 11,
    temperature_unit=celsius,
    times=times.tolist(),
    time_unit=minute,
)

from devtools import pprint

pprint(p)
