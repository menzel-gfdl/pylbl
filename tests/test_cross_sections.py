from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.hitran.cross_sections import HitranCrossSection


class TestCrossSections(TestCase):

    def test_cross_sections(self):
        formulae = ["CFC-11", "CFC-12"]
        grid = arange(1., 3250., 0.01)
        for formula in formulae:
            molecule = HitranCrossSection(formula)
            for pressure in [10000., 100000.]:
                for temperature in [211.11, 250.23, 290.11]:
                    molecule.absorption_coefficient(temperature, pressure, grid)


if __name__ == "__main__":
    main()
