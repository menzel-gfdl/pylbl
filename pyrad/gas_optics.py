from .lbl.hitran.lines import SpectralLines
from .lbl.tips import TotalPartitionFunction


class Gas(object):
    def __init__(self, formula, database, isotopologues="all"):
        partition_function = TotalPartitionFunction()
        partition_function.download_from_web(formula)
        self.spectral_lines = SpectralLines(formula, isotopologues, database,
                                            partition_function)

    def absorption_coefficient(self, temperature, pressure, partial_pressure, spectral_grid,
                               line_profile, cut_off=25.):
        return self.spectral_lines.absorption_coefficient(temperature, pressure, partial_pressure,
                                                          spectral_grid, line_profile, cut_off)
