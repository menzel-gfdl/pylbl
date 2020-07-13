from logging import basicConfig, INFO
from unittest import main, TestCase

from numpy import seterr

from pyrad.lbl.hitran.database import Database
from pyrad.lbl.hitran.lines import SpectralLines
from pyrad.lbl.line_profiles import Voigt
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.utils.grids import UniformGrid1D


class TestLines(TestCase):

    def test_lines(self):
        #Configure layer.
        mb_to_atm = 0.000986923
        spectral_grid = UniformGrid1D(1., 3000., 0.1)

        #Surface layer of CIRC case 1.
        temperature = 288.99 #K
        pressure = 983.88*mb_to_atm #atm
        abundance = {"H2O" : 0.006637074}
        hitran = Database()
        profile = Voigt()
        tips = TotalPartitionFunction()
        for formula, concentration in abundance.items():
            partial_pressure = pressure*concentration #atm
            tips.download_from_web(formula)
            lines = SpectralLines(formula, "all", hitran, tips)
            k = lines.absorption_coefficient(temperature, pressure, partial_pressure,
                                             spectral_grid.points, profile)


if __name__ == "__main__":
    basicConfig(format="%(asctime)-15s - %(pathname)s(%(lineno)d):\n\t%(message)s",
                level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    main()
