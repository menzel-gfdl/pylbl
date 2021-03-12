from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.hitran.cross_sections import HitranCrossSection


formulae = ["CFC-11", "CFC-12"]
grid = arange(1., 3250., 0.01)
pressures = [10000., 100000.]
temperatures = [211.11, 250.23, 290.11]


class TestCrossSections(TestCase):

    def test_cross_sections(self):
        for formula in formulae:
            molecule = HitranCrossSection(formula)
            for pressure in pressures:
                for temperature in temperatures:
                    molecule.absorption_coefficient(temperature, pressure, grid)

    def test_cross_sections_database(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            for formula in formulae:
                molecule = HitranCrossSection(formula)
                molecule.create_database(db)
            for formula in formulae:
                molecule = HitranCrossSection(formula, database=db)
                for pressure in pressures:
                    for temperature in temperatures:
                        molecule.absorption_coefficient(temperature, pressure, grid)

    def test_results(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            for formula in formulae:
                molecule_net = HitranCrossSection(formula)
                molecule_net.create_database(db)
                molecule_db = HitranCrossSection(formula, database=db)
                for band_net, band_db in zip(molecule_net.bands, molecule_db.bands):
                    for table_net, table_db in zip(band_net, band_db):
                        if table_net.band_params != table_db.band_params or \
                           table_net.temperature != table_db.temperature or \
                           table_net.pressure != table_db.pressure:
                            raise ValueError("Non-identical tables (download vs database).")
                for pressure in pressures:
                    for temperature in temperatures:
                        x_net = molecule_net.absorption_coefficient(temperature, pressure, grid)
                        x_db = molecule_db.absorption_coefficient(temperature, pressure, grid)
                        if x_net.tobytes() != x_db.tobytes():
                            raise ValueError("Non-identical cross sections (downloaded vs database).")


if __name__ == "__main__":
    main()
