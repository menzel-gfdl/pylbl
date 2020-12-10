from os.path import join
from unittest import main, TestCase

from pyrad.optics.clouds.liquid import LiquidCloudOptics
from pyrad.utils.grids import UniformGrid1D


class TestTips(TestCase):

    def test_hu_stamnes(self):
        clouds = LiquidCloudOptics(join("pyrad_data", "clouds", "hu_stamnes.nc"))
        clouds.optics(water_content=50., equivalent_radius=10.,
                      grid=UniformGrid1D(5000., 6000., 10))


if __name__ == "__main__":
    main()
