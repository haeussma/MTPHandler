from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate


class MTPReaderFactory:
    @staticmethod
    def read(
        cls: Plate,
        path: str,
        ph: Optional[float] = None,
        wavelength: Optional[float] = None,
        time: Optional[List[float]] = None,
        time_unit: Optional[str] = None,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
    ) -> Plate:
        from MTPHandler.readers.megellan_parser import read_magellan
        from MTPHandler.readers.multiskan_spectrum_parser import read_multiskan_spectrum
        from MTPHandler.readers.spectramax_parser import read_spectramax
        from MTPHandler.readers.tekan_spark_parser import read_tekan_spark

        try:
            return read_magellan(
                cls=cls,
                path=path,
                ph=ph,
                wavelength=wavelength,
            )
        except Exception:
            pass

        try:
            return read_multiskan_spectrum(
                cls=cls,
                path=path,
                time=time,
                time_unit=time_unit,
                ph=ph,
                temperature=temperature,
                temperature_unit=temperature_unit,
            )
        except Exception:
            pass

        try:
            return read_spectramax(
                cls=cls,
                path=path,
                ph=ph,
                time_unit=time_unit,
            )
        except Exception:
            pass

        try:
            return read_tekan_spark(
                cls=cls,
                path=path,
                ph=ph,
            )
        except Exception:
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
