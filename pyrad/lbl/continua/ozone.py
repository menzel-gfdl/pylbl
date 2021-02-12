from .utils import cm_to_m, download_gfdl_data


class OzoneContinuum(object):
    """Ozone continuum.

    Attributes:
        cross_section: Cross section [cm2].
    """
    def __init__(self, database=None):
        if database is None:
            self.download()

    def download(self):
        """Downloads cross sections from GFDL's FTP site."""
        download_gfdl_data(self, "ozone_continuum", ["ozone_continuum.csv",],
                           ["cross_section",])

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates the absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Numpy array of wavenumbers [cm-1].

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        return self.cross_section.regrid(grid)[:]*cm_to_m*cm_to_m
