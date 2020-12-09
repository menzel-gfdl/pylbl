from logging import basicConfig, INFO
from os.path import join
from unittest import main, TestCase

from numpy import seterr

from pyrad.optics.aerosols import AerosolOptics
from pyrad.utils.grids import UniformGrid1D


class TestAerosols(TestCase):

    def test_aerosols(self):
        sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
        sulfate.optics(concentration=1., grid=UniformGrid1D(1., 3000., .1),
                       humidity=50.5, mixture=80.3)


if __name__ == "__main__":
    basicConfig(level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    main()
