from numpy import log, pi, sqrt
from scipy.special import wofz

from .doppler import doppler_broadened_halfwidth
from .line_parameters import PARAMETERS
from .lorentz import pressure_broadened_halfwidth


class Voigt(object):
    def __init__(self):
        self.parameters = [PARAMETERS[x] for x in ["center", "gamma_air", "gamma_self", "n_air"]]
        self.doppler_halfwidth = None
        self.pressure_halfwidth = None

    def update(self, spectral_lines, temperature, pressure, partial_pressure):
        """Calculate per-spectral-line pressure-broadened halfwidths.

        Args:
            spectral_lines: SpectralLines object.
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
        """
        self.doppler_halfwidth = doppler_broadened_halfwidth(temperature, spectral_lines.mass,
                                                             spectral_lines.v)
        self.pressure_halfwidth = pressure_broadened_halfwidth(pressure, partial_pressure, temperature,
                                                               spectral_lines.n_air, spectral_lines.gamma_air,
                                                               spectral_lines.gamma_self)

    def profile(self, spectral_lines, v, index):
        """Calculate Voigt profiles.

        Args:
            spectral_lines: SpectralLines object.
            v: Wavenumber [cm-1].
            index: Spectral line index.
        """
        return voigt_profile(v - spectral_lines.v[index], self.pressure_halfwidth[index],
                             self.doppler_halfwidth[index])


def voigt_profile(dv, pressure_halfwidth, doppler_halfwidth):
    """Calculates a Voigt line profile.

    Args:
        dv: Wavenumber distance from line center [cm-1].
        pressure_halfwidth: Pressure-broadened line half-width [cm -1].
        doppler_halfwidth: Doppler-broadened line half-width [cm -1].

    Returns:
        Voigt line profile [cm-1].
    """
    sigma = doppler_halfwidth/sqrt(2*log(2))
    return wofz((dv + 1j*pressure_halfwidth)/sigma/sqrt(2)).real/sigma/sqrt(2*pi)
