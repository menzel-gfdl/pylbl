from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from hitran import create_database, HitranDatabase, HitranSpectralLines, write_to_ascii
from tips import TotalPartitionFunction


name = "test"


class TestHitranHTMLParsers(TestCase):

    def test_download_from_web(self):
        molecule = "H2O"
        total_partition_function = TotalPartitionFunction()
        total_partition_function.download_from_web(molecule)
        spectral_database = HitranDatabase()
        parameters = spectral_database.download_from_web(molecule, upper_bound=1000.)
        isotopologues = spectral_database.isotopologues[molecule]
        spectral_lines = HitranSpectralLines(parameters, isotopologues, total_partition_function)

    def test_read_from_file(self):
        molecule = "CO2"
        total_partition_function = TotalPartitionFunction()
        total_partition_function.download_from_web(molecule)
        spectral_database = HitranDatabase()
        parameters = spectral_database.download_from_web(molecule, upper_bound=1000.)
        with TemporaryDirectory() as directory:
            id = spectral_database.molecule_ids[molecule]
            path = join(directory, name)
            write_to_ascii(path, id, parameters)
            parameters = spectral_database.read_from_ascii(path, molecule, upper_bound=1000.)
        isotopologues = spectral_database.isotopologues[molecule]
        spectral_lines = HitranSpectralLines(parameters, isotopologues, total_partition_function)

    def test_from_database(self):
        molecule = "H2O"
        spectral_database = HitranDatabase()
        parameters = spectral_database.download_from_web(molecule, upper_bound=1000.)
        with TemporaryDirectory() as directory:
            database = join(directory, "{}.db".format(name))
            create_database(database, molecule, parameters)
            spectral_database.load_from_database(molecule, database)


if __name__ == "__main__":
    main()
