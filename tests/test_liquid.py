from os.path import join
from unittest import main, TestCase

from pyrad.optics.clouds.liquid import LiquidCloudOptics
from pyrad.utils.grids import UniformGrid1D


class TestTips(TestCase):

    def test_hu_stamnes(self):
        clouds = LiquidCloudOptics(join("pyrad_data", "clouds", "hu_stamnes.nc"))
        lwc =  50. #g m-3
        r = 10. #micron
        grid = UniformGrid1D(5000., 6000., 10) #cm-1
        optics = clouds.optics(lwc, r, grid)


if __name__ == "__main__":
    main()
