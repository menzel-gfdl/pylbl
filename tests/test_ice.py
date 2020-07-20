from os.path import join
from unittest import main, TestCase

from pyrad.optics.clouds.ice import IceCloudOptics
from pyrad.utils.grids import UniformGrid1D


class TestTips(TestCase):

    def test_fu_liou_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "fu_liou.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="longwave")

    def test_ebert_curry_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "ebert_curry.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="longwave")

    def test_chou_suarez_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="longwave")

    def test_fu_liou_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "fu_liou.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="shortwave")

    def test_ebert_curry_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "ebert_curry.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="shortwave")

    def test_chou_suarez_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
        iwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(500., 1500., 10)
        optics = clouds.optics(iwc, r, grid, mode="shortwave")


if __name__ == "__main__":
    main()
