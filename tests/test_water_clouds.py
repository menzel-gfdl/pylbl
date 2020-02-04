from os.path import join
from unittest import main, TestCase

from numpy import linspace

from pyrad.clouds.water import HuStamnes, Slingo


class TestTips(TestCase):

    def test_hu_stamnes(self):
        clouds = HuStamnes(join("pyrad_data", "clouds", "hu_stamnes.nc"))
        lwc =  50. #g m-3
        r = 10. #micron
        grid = linspace(5000., 6000., 10) #cm-1
        beta, omega, g = clouds.optics(lwc, r, grid)
        z = 1.e3 #cm

    def test_slingo(self):
        clouds = Slingo(join("pyrad_data", "clouds", "slingo.nc"))
        lwc =  50. #g m-3
        r = 10. #micron
        grid = linspace(5000., 6000., 10) #cm-1
        z = 1.e3 #cm
        cm_to_m = 1.e-2 #[m cm-1]
        lwp = lwc*z*cm_to_m
        tau, omega, g = clouds.optics(lwp, r, grid)


if __name__ == "__main__":
    main()
