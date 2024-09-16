from MTPHandler.dataclasses import Plate
from MTPHandler.readers import read_multiskan_spectrum
from MTPHandler.units import nm, s

def test_spectramax():

    path = "tests/data/ABTS_EnzymeML_340nm_420nm_2.5x_pH3_25deg.txt"
    ph = 6.9
    time = [0,1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15]
