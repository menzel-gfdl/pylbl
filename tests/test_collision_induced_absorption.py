from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA


grid = arange(1., 3250., 0.01)
#interactions = [["CO2", "CO2"], ["N2", "N2"], ["O2", "N2"], ["O2", "O2"]]
interactions = [["O2", "N2"],]
temperatures = [211.11, 235.45, 250.23, 270.45, 290.11, 300.1]


class TestCollisionInducedAbsorption(TestCase):

    def test_collision_induced_absorption(self):
        for interaction in interactions:
            molecule = HitranCIA(*interaction)
            for temperature in temperatures:
                molecule.absorption_coefficient(temperature, grid)

    def test_collision_induced_absorption_database(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            for interaction in interactions:
                molecule = HitranCIA(*interaction)
                molecule.create_database(db)
            for interaction in interactions:
                molecule = HitranCIA(*interaction, database=db)
                for temperature in temperatures:
                    molecule.absorption_coefficient(temperature, grid)

    def test_results(self):
        with TemporaryDirectory() as directory:
            db = join(directory, "test.sqlite")
            for interaction in interactions:
                molecule_net = HitranCIA(*interaction)
                molecule_net.create_database(db)
                molecule_db = HitranCIA(*interaction, database=db)
                for band_net, band_db in zip(molecule_net.bands, molecule_db.bands):
                    for table_net, table_db in zip(band_net, band_db):
                        if table_net.band_params != table_db.band_params or \
                           table_net.temperature != table_db.temperature:
                            raise ValueError("Non-identical tables (download vs database).")
                for temperature in temperatures:
                    x_net = molecule_net.absorption_coefficient(temperature, grid)
                    x_db = molecule_db.absorption_coefficient(temperature, grid)
                    if x_net.tobytes() != x_db.tobytes():
                        raise ValueError("Non-identical cross sections (downloaded vs database).")


if __name__ == "__main__":
    main()
