from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.continua import WaterVaporForeignContinuum, WaterVaporSelfContinuum


grid = arange(1., 3250., 0.01)
pressure = 101300.
temperature = 300.
volume_mixing_ratio = 0.02


class TestWaterVaporContinuum(TestCase):

    def test_foreign_continuum(self):
        WaterVaporForeignContinuum().absorption_coefficient(temperature, pressure,
                                                            volume_mixing_ratio, grid)

    def test_foreign_continuum_database(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            h2o = WaterVaporForeignContinuum()
            h2o.create_database(db)
            h2o = WaterVaporForeignContinuum(db)
            h2o.absorption_coefficient(temperature, pressure, volume_mixing_ratio, grid)

    def test_foreign_results(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            h2o_net = WaterVaporForeignContinuum()
            h2o_net.create_database(db)
            h2o_db = WaterVaporForeignContinuum(database=db)
            x_net = h2o_net.absorption_coefficient(temperature, pressure,
                                                   volume_mixing_ratio, grid)
            x_db = h2o_db.absorption_coefficient(temperature, pressure,
                                                 volume_mixing_ratio, grid)
            if x_net.tobytes() != x_db.tobytes():
                raise ValueError("Non-identical cross sections (downloaded vs database).")

    def test_self_continuum(self):
        WaterVaporSelfContinuum().absorption_coefficient(temperature, pressure,
                                                         volume_mixing_ratio, grid)

    def test_self_continuum_database(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            h2o = WaterVaporSelfContinuum()
            h2o.create_database(db)
            h2o = WaterVaporSelfContinuum(db)
            h2o.absorption_coefficient(temperature, pressure, volume_mixing_ratio, grid)

    def test_self_results(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            h2o_net = WaterVaporSelfContinuum()
            h2o_net.create_database(db)
            h2o_db = WaterVaporSelfContinuum(database=db)
            x_net = h2o_net.absorption_coefficient(temperature, pressure,
                                                   volume_mixing_ratio, grid)
            x_db = h2o_db.absorption_coefficient(temperature, pressure,
                                                 volume_mixing_ratio, grid)
            if x_net.tobytes() != x_db.tobytes():
                raise ValueError("Non-identical cross sections (downloaded vs database).")


if __name__ == "__main__":
    main()
