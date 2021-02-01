from numpy import log, pi, sqrt
from scipy.special import wofz

from .doppler import doppler_broadened_halfwidth
from .line_mixing import LineMixing
from .lorentz import pressure_broadened_halfwidth


class Voigt(object):
    """Voigt line profile.

    Attributes:
        doppler_halfwidth: Doppler-broadened halfwidth [cm-1].
        line_mixing: LineMixing object.
        parameters: List of HITRAN parameter names.
        pressure_halfwidth: Pressure-broadened halfwidth [cm-1].
    """

    def __init__(self, relaxation_matrix_path=None):
        self.parameters = ["center", "gamma_air", "gamma_self", "n_air", "n_self"]
        self.doppler_halfwidth = None
        self.pressure_halfwidth = None
        if relaxation_matrix_path is not None:
            self.parameters += ["upper_level_quantum_numbers",
                                "lower_level_quantum_numbers",
                                "lower_level_degeneracy"]
            self.line_mixing = LineMixing(relaxation_matrix_path)
        else:
            self.line_mixing = None

    def update(self, spectral_lines, temperature, pressure, partial_pressure):
        """Calculate per-spectral-line pressure-broadened halfwidths.

        Args:
            spectral_lines: SpectralLines object.
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
        """
        if self.line_mixing is not None:
            self.y = self.line_mixing.first_order_coefficients(spectral_lines,
                                                               temperature,
                                                               pressure,
                                                               partial_pressure/pressure)
        self.doppler_halfwidth = doppler_broadened_halfwidth(temperature, spectral_lines.mass,
                                                             spectral_lines.v)
        self.pressure_halfwidth = pressure_broadened_halfwidth(pressure, partial_pressure,
                                                               temperature,
                                                               spectral_lines.n_air,
                                                               spectral_lines.n_self,
                                                               spectral_lines.gamma_air,
                                                               spectral_lines.gamma_self)

    def profile(self, spectral_lines, v, index):
        """Calculate Voigt profiles.

        Args:
            spectral_lines: SpectralLines object.
            v: Wavenumber [cm-1].
            index: Spectral line index.

        Returns:
            Line broadening [cm].
        """
        x = voigt_profile(v - spectral_lines.v[index], self.pressure_halfwidth[index],
                          self.doppler_halfwidth[index])
        if self.line_mixing is None:
            return x.real
        else:
            return ((1. - 1j*self.y[index])*x).real


def voigt_profile(dv, pressure_halfwidth, doppler_halfwidth):
    """Calculates a Voigt line profile.

    Args:
        dv: Wavenumber distance from line center [cm-1].
        pressure_halfwidth: Pressure-broadened line half-width [cm -1].
        doppler_halfwidth: Doppler-broadened line half-width [cm -1].

    Returns:
        Voigt line profile broadening [cm].
    """
    sigma = doppler_halfwidth/sqrt(2*log(2))
    return wofz((dv + 1j*pressure_halfwidth)/sigma/sqrt(2))/sigma/sqrt(2*pi)
