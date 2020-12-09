from os.path import join
from unittest import main, TestCase

from pyrad.optics.clouds.ice import IceCloudOptics
from pyrad.utils.grids import UniformGrid1D


class TestTips(TestCase):

    def test_fu_liou_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "fu_liou.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="longwave")

    def test_ebert_curry_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "ebert_curry.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="longwave")

    def test_chou_suarez_lw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="longwave")

    def test_fu_liou_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "fu_liou.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="shortwave")

    def test_ebert_curry_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "ebert_curry.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="shortwave")

    def test_chou_suarez_sw(self):
        clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
        clouds.optics(ice_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(500., 1500., 10), mode="shortwave")


if __name__ == "__main__":
    main()
