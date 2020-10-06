from logging import basicConfig, INFO
from unittest import main, TestCase

from numpy import seterr

from pyrad.optics.gas import Gas
from pyrad.utils.grids import UniformGrid1D


class TestGasOptics(TestCase):

    def test_gas_optics(self):
        #Configure layer.
        mb_to_atm = 0.000986923
        spectral_grid = UniformGrid1D(1., 3000., 0.1)

        #Surface layer of CIRC case 1.
        temperature = 288.99 #K
        pressure = 983.88*mb_to_atm #atm
        abundance = {"H2O" : 0.006637074}
        for formula, concentration in abundance.items():
            partial_pressure = pressure*concentration #atm
            gas = Gas(formula)
            k = gas.absorption_coefficient(temperature, pressure, partial_pressure,
                                           spectral_grid.points)


if __name__ == "__main__":
    basicConfig(format="%(asctime)-15s - %(pathname)s(%(lineno)d):\n\t%(message)s",
                level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    main()
