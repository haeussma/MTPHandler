import numpy as np

from mtphandler.model import Plate
from mtphandler.units import celsius, second

times = np.linspace(0, 10, 11)


p = Plate(
    id="MTP_001",
    name="Enzyme Kinetics",
    temperatures=[25] * 11,
    temperature_unit=celsius,
    times=times.tolist(),
    time_unit=second,
)

from devtools import pprint

pprint(p)
