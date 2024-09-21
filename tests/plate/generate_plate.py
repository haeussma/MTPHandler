import numpy as np

<<<<<<< HEAD
from mtphandler.model import Plate
from mtphandler.units import celsius, second
=======
from MTPHandler.model import Plate
from MTPHandler.units import celsius, minute
>>>>>>> main

times = np.linspace(0, 10, 11)


p = Plate(
    id="MTP_001",
    name="Enzyme Kinetics",
    temperatures=[25] * 11,
    temperature_unit=celsius,
    times=times.tolist(),
<<<<<<< HEAD
    time_unit=second,
=======
    time_unit=minute,
>>>>>>> main
)

from devtools import pprint

pprint(p)
