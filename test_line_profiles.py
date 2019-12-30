from unittest import main, TestCase

from numpy import linspace, zeros

from line_profiles import Doppler, Lorentz, Voigt


grid = linspace(-1., 1., 1000)
halfwidth = 0.1
v_center = 0.
tolerance = 1.e-8

class TestLineProfiles(TestCase):

    def test_doppler(self):
        profile = Doppler()
        profile.alpha = halfwidth
        k = zeros(grid.size)
        for i, v in enumerate(grid):
            k[i] = profile.profile(v, v_center)

    def test_lorentz(self):
        profile = Lorentz()
        profile.gamma = halfwidth
        k = zeros(grid.size)
        for i, v in enumerate(grid):
            k[i] = profile.profile(v, v_center)

    def test_voigt(self):
        profile = Voigt()
        profile.alpha = halfwidth
        profile.gamma = halfwidth
        k = zeros(grid.size)
        for i, v in enumerate(grid):
            k[i] = profile.profile(v, v_center)


if __name__ == "__main__":
    main()
