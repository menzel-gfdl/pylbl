from ..lbl.hitran import Hitran, Voigt
from ..lbl.tips import TotalPartitionFunction


class Gas(object):
    def __init__(self, formula, hitran_database=None, isotopologues=None, 
                 line_profile=Voigt(), tips_database=None):
        database = Hitran(formula, line_profile, isotopologues, hitran_database)
        partition_function = TotalPartitionFunction(formula, tips_database)
        self.spectral_lines = database.spectral_lines(partition_function)

    def absorption_coefficient(self, temperature, pressure, partial_pressure, spectral_grid,
                               line_cut_off=25.):
        """Calculates absorption coefficients for the gas using line-by-line method.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial Pressure [atm].
            spectral_grid: Wavenumber grid [cm-1].
            line_cut_off: Cut-off from spectral line center [cm-1].

        Returns:
            Absorption coefficients [cm-1].
        """
        return self.spectral_lines.absorption_coefficient(temperature, pressure, partial_pressure,
                                                          spectral_grid, line_cut_off)
