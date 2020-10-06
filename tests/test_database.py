from unittest import main, TestCase

from pyrad.lbl.hitran import Hitran, Doppler, Lorentz, Voigt


class TestDatabase(TestCase):

    def test_database_doppler(self):
        formulae = ["H2O", "CO2", "O3", "N2O", "CH4"]
        for formula in formulae:
            hitran = Hitran(formula, Doppler())

    def test_database_lorentz(self):
        formulae = ["H2O", "CO2", "O3", "N2O", "CH4"]
        for formula in formulae:
            hitran = Hitran(formula, Lorentz())

    def test_database_voigt(self):
        formulae = ["H2O", "CO2", "O3", "N2O", "CH4"]
        for formula in formulae:
            hitran = Hitran(formula, Voigt())


if __name__ == "__main__":
    main()
