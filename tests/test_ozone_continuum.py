from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.continua import OzoneContinuum


class TestOzoneContinuum(TestCase):

    def test_continuum(self):
        OzoneContinuum().absorption_coefficient(arange(9000., 50000., 0.1))

if __name__ == "__main__":
    main()
