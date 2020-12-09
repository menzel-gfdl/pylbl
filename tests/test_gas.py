from logging import basicConfig, INFO
from unittest import main, TestCase

from numpy import seterr

from pyrad.optics.gas import Gas
from pyrad.utils.grids import UniformGrid1D


class TestGasOptics(TestCase):

    def test_gas_optics(self):
        # Surface layer of CIRC case 1.
        abundance = {"H2O": 0.006637074}
        for formula, concentration in abundance.items():
            gas = Gas(formula)
            gas.absorption_coefficient(temperature=288.99, pressure=98388.,
                                       volume_mixing_ratio=concentration,
                                       spectral_grid=UniformGrid1D(1., 3000., 0.1).points)


if __name__ == "__main__":
    basicConfig(format="%(asctime)-15s - %(pathname)s(%(lineno)d):\n\t%(message)s",
                level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    main()
