from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.continua import WaterVaporForeignContinuum, WaterVaporSelfContinuum


class TestWaterVaporContinuum(TestCase):

    def test_foreign_continuum(self):
        WaterVaporForeignContinuum().absorption_coefficient(300., 101300., 0.02,
                                                            arange(1., 3250., 0.01))
    def test_self_continuum(self):
        WaterVaporSelfContinuum().absorption_coefficient(300., 101300., 0.02,
                                                         arange(1., 3250., 0.01))

if __name__ == "__main__":
    main()
