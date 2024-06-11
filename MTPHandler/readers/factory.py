from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from MTPHandler.core.plate import Plate


class MTPReaderFactory:
    @staticmethod
    def read(
        cls: Plate,
        path: str,
        ph: float|None = None,
        wavelength: float|None = None,
        time: list[float]|None = None,
        time_unit: str|None = None,
        temperature: float|None = None,
        temperature_unit: str|None = None,
    ) -> Plate:
        from MTPHandler.readers.megellan_parser import read_magellan
        from MTPHandler.readers.multiskan_spectrum_parser import read_multiskan_spectrum
        from MTPHandler.readers.spectramax_parser import read_spectramax
        from MTPHandler.readers.tekan_spark_parser import read_tekan_spark
        from MTPHandler.readers.biotek import read_biotek


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

        try:
            return read_biotek(
                cls=cls,
                path=path,
                ph=ph,
            )
        except Exception:
            pass

        raise ValueError("Could not read file with implemented readers.")
