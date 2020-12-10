from numpy import exp, log, pi, sqrt


class Doppler(object):
    """Doppler line profile.

    Attributes:
        halfwidth: Doppler-broadened halfwidth [cm-1].
        parameters: List of HITRAN parameter names.
    """

    def __init__(self):
        self.parameters = ["center"]
        self.halfwidth = None

    def update(self, spectral_lines, temperature, *args, **kwargs):
        """Calculate per-spectral-line doppler-broadened halfwidths.

        Args:
            spectral_lines: SpectralLines object.
            temperature: Temperature [K].
        """
        self.halfwidth = doppler_broadened_halfwidth(temperature, spectral_lines.mass,
                                                     spectral_lines.v)

    def profile(self, spectral_lines, v, index):
        """Calculate doppler profiles.

        Args:
            spectral_lines: SpectralLines object.
            v: Wavenumber [cm-1].
            index: Spectral line index.

        Returns:
            Line broadening [cm].
        """
        return doppler_profile(v - spectral_lines.v[index], self.halfwidth[index])


def doppler_broadened_halfwidth(temperature, mass, transition_wavenumber):
    """Calculate the doppler-broadened line halfwidth.

    Args:
        temperature: Temperature [K].
        mass: Molecular mass [g].
        transition_wavenumber: Transition wavenumber [cm-1].

    Returns:
        Doppler-broadened line halfwidth [cm-1].
    """
    m = mass/6.023e23  # Mass/Avagadro's number [g].
    c = 2.99792458e10  # Speed of light [cm s-1].
    kb = 1.380658e-16  # Boltzmann constant [erg K-1].
    return sqrt(log(2.))*transition_wavenumber*sqrt(2.*kb*temperature/(m*c*c))


def doppler_profile(dv, halfwidth):
    """Calculates a Doppler line profile.

    Args:
        dv: Wavenumber distance from line center [cm-1].
        halfwidth: Doppler-broadened line half-width [cm -1].

    Returns:
        Doppler line profile broadening [cm].
    """
    alpha_d = halfwidth/sqrt(log(2))
    return exp(-1.*dv*dv/(alpha_d*alpha_d))/(alpha_d*sqrt(pi))
