from unittest import main, TestCase

from pyrad.lbl.hitran import Hitran, Doppler, Lorentz, Voigt


class TestDatabase(TestCase):

    def test_database_doppler(self):
        formulae = ["H2O"]
        for formula in formulae:
            Hitran(formula, Doppler())

    def test_database_lorentz(self):
        formulae = ["H2O"]
        for formula in formulae:
            Hitran(formula, Lorentz())

    def test_database_voigt(self):
        formulae = ["H2O"]
        for formula in formulae:
            Hitran(formula, Voigt())


if __name__ == "__main__":
    main()
