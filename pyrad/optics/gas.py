from ..lbl.continua import OzoneContinuum, WaterVaporContinuum
from ..lbl.hitran import Hitran, Voigt
from ..lbl.tips import TotalPartitionFunction


pa_to_atm = 9.86923e-6  # [atm Pa-1].
cm_to_m = 0.01  # [m cm-1].


class Gas(object):
    """Line-by-line gas optics.

    Attributes:
        spectral_lines: SpectralLines object.
    """
    def __init__(self, formula, hitran_database=None, isotopologues=None,
                 line_profile=Voigt(), tips_database=None):
        """Obtains molecular line parameters.

        Args:
            formula: Chemical formula.
            hitran_database: Path to sqlite hitran database.
            isotopologues: Lists of Isotopologue objects.
            line_profile: Doppler, Lorentz, or Voigt object.
            tips_database: Path to sqlite tips database.
        """
        database = Hitran(formula, line_profile, isotopologues, hitran_database)
        partition_function = TotalPartitionFunction(formula, tips_database)
        self.spectral_lines = database.spectral_lines(partition_function)
        if formula == "H2O":
            self.continuum = WaterVaporContinuum()
        elif formula == "O3":
            self.continuum = OzoneContinuum()

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio,
                               spectral_grid, line_cut_off=25.):
        """Calculates absorption coefficients for the gas using line-by-line method.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            spectral_grid: Wavenumber grid [cm-1].
            line_cut_off: Cut-off from spectral line center [cm-1].

        Returns:
            Absorption coefficients [m2].
        """
        lines = self.spectral_lines.absorption_coefficient(temperature, pressure*pa_to_atm,
                                                           pressure*pa_to_atm*volume_mixing_ratio,
                                                           spectral_grid, line_cut_off) * \
            cm_to_m*cm_to_m
        if "continuum" in vars(self):
            continuum = self.continuum.absorption_coefficient(temperature, pressure,
                                                              volume_mixing_ratio, spectral_grid)
            return lines + continuum
        else:
            return lines
