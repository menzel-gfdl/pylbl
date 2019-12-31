from tempfile import NamedTemporaryFile
from unittest import main, TestCase

from pylbl.hitran import create_database, HitranDatabase, HitranSpectralLines, write_to_ascii
from pylbl.tips import TotalPartitionFunction


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
        with NamedTemporaryFile() as path:
            id = spectral_database.molecule_ids[molecule]
            write_to_ascii(path.name, id, parameters)
            parameters = spectral_database.read_from_ascii(path.name, molecule, upper_bound=1000.)
        isotopologues = spectral_database.isotopologues[molecule]
        spectral_lines = HitranSpectralLines(parameters, isotopologues, total_partition_function)

    def test_from_database(self):
        molecule = "H2O"
        spectral_database = HitranDatabase()
        parameters = spectral_database.download_from_web(molecule, upper_bound=1000.)
        with NamedTemporaryFile(suffix=".db") as database:
            create_database(database.name, molecule, parameters)
            spectral_database.load_from_database(molecule, database.name)


if __name__ == "__main__":
    main()
