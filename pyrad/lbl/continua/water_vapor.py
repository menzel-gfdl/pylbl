from numpy import exp

from .utils import cm_to_m, download_gfdl_data, Pa_to_atm

class _WaterVaporContinuum(object):
    """Water vapor self continuum.

    Attributes:
        c: Broadening coefficient [cm2 atm-1].
        t: Broadening temperature coefficient [K-1].
    """
    def __init__(self, database=None):
        if database is None:
            self.download()
        else:
            raise NotImplementedError("database reads not implemented yet.")

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
        partial_pressure = self.partial_pressure(pressure, volume_mixing_ratio)
        return (296./temperature)*self.c.regrid(grid)[:]*partial_pressure * \
               exp(self.t.regrid(grid)[:]*(296. - temperature))*cm_to_m*cm_to_m

    def download(self, filenames):
        """Downloads water vapor continuum coefficients from GFDL's FTP site."""
        parameters = ["c_self", "t_self"]
        download_gfdl_data(self, "water_vapor_continuum", filenames, ["c", "t"])

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        raise NotImplementedError("you must override this function.")


class WaterVaporSelfContinuum(_WaterVaporContinuum):
    """Water vapor self continuum."""
    def download(self):
        """Downloads water vapor self-continuum coefficients from GFDL's FTP site."""
        super().download(["296MTCKD25_S.csv", "CKDS.csv"])

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        return pressure*volume_mixing_ratio*Pa_to_atm


class WaterVaporForeignContinuum(_WaterVaporContinuum):
    """Water vapor foreign continuum."""
    def download(self):
        """Downloads water vapor continuum foreign-coefficients from GFDL's FTP site."""
        super().download(["296MTCKD25_F.csv", "CKDF.csv"])

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        return pressure*(1. - volume_mixing_ratio)*Pa_to_atm
