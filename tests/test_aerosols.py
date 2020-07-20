from logging import basicConfig, INFO
from os.path import join
from unittest import main, TestCase

from numpy import seterr

from pyrad.optics.aerosols import AerosolOptics
from pyrad.utils.grids import UniformGrid1D


class TestAerosols(TestCase):

    def test_aerosols(self):
        grid = UniformGrid1D(1., 3000., .1)
        concentration = 1. #
        relative_humidity = 50.5 #[%].
        mixture = 80.3 #[%].
        sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
        optics = sulfate.optics(concentration, grid, relative_humidity, mixture)


if __name__ == "__main__":
    basicConfig(level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    main()
