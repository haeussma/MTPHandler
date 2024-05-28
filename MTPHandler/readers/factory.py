from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate


class MTPReaderFactory:
    @staticmethod
    def read(
        path: str,
        ph: Optional[float] = None,
        wavelength: Optional[float] = None,
        time: Optional[List[float]] = None,
        time_unit: Optional[str] = None,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
    ) -> Plate:
        from MTPHandler.core import Plate
        from MTPHandler.readers.megellan_parser import read_magellan
        from MTPHandler.readers.multiskan_parser import read_multiskan
        from MTPHandler.readers.spectramax_parser import read_spectramax

        try:
            return read_magellan(Plate, path, ph, wavelength)
        except Exception:
            pass

        try:
            return read_multiskan(
                Plate, path, time, time_unit, ph, temperature, temperature_unit
            )
        except ValueError:
            pass

        try:
            return read_spectramax(Plate, path, ph, time_unit)
        except ValueError:
            pass

        raise ValueError("Could not read file with implemented readers.")


if __name__ == "__main__":
    path = "tests/data/magellan.xlsx"

    ph = 7.0
    wavelength = 450.0

    plate = MTPReaderFactory.read(
        path=path,
        ph=ph,
        wavelength=wavelength,
    )

    print(plate)
