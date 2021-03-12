from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.continua import OzoneContinuum


grid = arange(9000., 50000., 0.1)


class TestOzoneContinuum(TestCase):

    def test_continuum(self):
        OzoneContinuum().absorption_coefficient(grid)

    def test_continuum_database(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            o3 = OzoneContinuum()
            o3.create_database(db)
            o3 = OzoneContinuum(database=db)
            o3.absorption_coefficient(grid)

    def test_results(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            o3_net = OzoneContinuum()
            o3_net.create_database(db)
            o3_db = OzoneContinuum(database=db)
            x_net = o3_net.absorption_coefficient(grid)
            x_db = o3_db.absorption_coefficient(grid)
            if x_net.tobytes() != x_db.tobytes():
                raise ValueError("Non-identical cross sections (downloaded vs database).")


if __name__ == "__main__":
    main()
