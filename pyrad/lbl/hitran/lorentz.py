from numpy import pi, power


class Lorentz(object):
    """Lorentz line profile.

    Attributes:
        halfwidth: Pressure-broadened halfwidth [cm-1].
        parameters: List of HITRAN parameter names.
    """

    def __init__(self):
        self.parameters = ["gamma_air", "gamma_self", "n_air"]
        self.halfwidth = None

    def update(self, spectral_lines, temperature, pressure, partial_pressure):
        """Calculate per-spectral-line pressure-broadened halfwidths.

        Args:
            spectral_lines: SpectralLines object.
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
        """
        self.halfwidth = pressure_broadened_halfwidth(pressure, partial_pressure, temperature,
                                                      spectral_lines.n_air, spectral_lines.gamma_air,
                                                      spectral_lines.gamma_self)

    def profile(self, spectral_lines, v, index):
        """Calculate lorentz profiles.

        Args:
            spectral_lines: SpectralLines object.
            v: Wavenumber [cm-1].
            index: Spectral line index.

        Returns:
            Line broadening [cm].
        """
        return lorentz_profile(v - spectral_lines.v[index], self.halfwidth[index])


def pressure_broadened_halfwidth(pressure, partial_pressure, temperature,
                                 n, gamma_air, gamma_self):
    """Calculates pressure-broadened line halfwidth.

    Args:
        pressure: Pressure [atm].
        partial_pressure: Partial pressure [atm].
        temperature: Temperature [K]
        n: Air-broadened temperature dependence powers.
        gamma_air: Air-broadened halfwidth [cm-1 atm-1].
        gamma_self: Self-broadened halfwidth [cm-1 atm-1].

    Returns:
        Pressure-broadened line halfwidth [cm-1].
    """
    return power((296./temperature), n) * \
        (gamma_air*(pressure - partial_pressure) + gamma_self*partial_pressure)


def lorentz_profile(dv, halfwidth):
    """Calculates a Lorentzian line profile.

    Args:
        dv: Wavenumber distance from line center [cm-1].
        halfwidth: Pressure-broadened line half-width [cm -1].

    Returns:
        Lorentz line profile broadening [cm].
    """
    return halfwidth/(pi*(halfwidth*halfwidth + dv*dv))
