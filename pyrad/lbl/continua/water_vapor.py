from numpy import exp

from .utils import cm_to_m, download_gfdl_data, Pa_to_atm


class WaterVaporContinuum(object):
    """Water vapor continuum.

    Attributes:
        c_foreign: Foreign-broadening coefficient [cm2 atm-1].
        c_self: Self-broadening coefficient [cm2 atm-1].
        t_foreign: Foreign-broadening temperature coefficient [K-1].
        t_foreign: Self-broadening temperature coefficient [K-1].
    """
    def __init__(self, database=None):
        if database is None:
            self.download()

    def download(self):
        """Downloads water vapor continuum coefficients from GFDL's FTP site."""
        names = ["296MTCKD25_F.csv", "296MTCKD25_S.csv", "CKDF.csv", "CKDS.csv"]
        parameters = ["c_foreign", "c_self", "t_foreign", "t_self"]
        download_gfdl_data(self, "water_vapor_continuum", names, parameters)

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
        pressure *= Pa_to_atm
        partial_pressure = pressure*volume_mixing_ratio
        k = (296./temperature)*((self.c_self.regrid(grid)[:]*partial_pressure *
                                 exp(self.t_self.regrid(grid)[:]*(296. - temperature))) +
                                (self.c_foreign.regrid(grid)[:]*(pressure - partial_pressure) *
                                 exp(self.t_foreign.regrid(grid)[:]*(296. - temperature))))
        return k*cm_to_m*cm_to_m
