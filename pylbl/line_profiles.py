from numpy import exp, log, pi, real, sqrt
from scipy.special import wofz


class LineProfile(object):
    """Base class for line profile objects.

    Attributes:
        gamma: Pressure-broadening (Lorentzian) halfwidth.
        alpha: Doppler-broadening (Gaussian) halfwidth.
    """
    def __init__(self):
        self.gamma = None
        self.alpha = None

    def profile(self, v, v_center):
        """Override this function with one that calculates the line profile.

        Args:
            v: Wavenumber.
            v_center: Transition wavenumber.

        Returns:
            Line profile.
        """
        raise NotImplementedError("not yet implemented.")


class Doppler(LineProfile):
    def profile(self, v, v_center):
        """Calculates a Doppler line profile.

        Args:
            v: Wavenumber [cm-1].
            v_center: Transition wavenumber [cm-1].

        Returns:
            Doppler line profile [cm-1].
        """
        alpha_d = self.alpha/sqrt(log(2))
        return exp(-1.*(v - v_center)*(v - v_center)/(alpha_d*alpha_d))/(alpha_d*sqrt(pi))


class Lorentz(LineProfile):
    def profile(self, v, v_center):
        """Calculates a Lorentzian line profile.

        Args:
            v: Wavenumber [cm-1].
            v_center: Transition wavenumber [cm-1].

        Returns:
            Lorentz line profile [cm-1].
        """
        return self.gamma/(pi*(self.gamma*self.gamma + (v - v_center)*(v - v_center)))


class Voigt(LineProfile):
    def profile(self, v, v_center):
        """Calculates a Voigt line profile.

        Args:
            v: Wavenumber [cm-1].
            v_center: Transition wavenumber [cm-1].

        Returns:
            Voigt line profile [cm-1].
        """
        sigma = self.alpha/sqrt(2*log(2))
        return real(wofz((v - v_center + 1j*self.gamma)/sigma/sqrt(2)))/sigma/sqrt(2*pi)
