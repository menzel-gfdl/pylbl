from logging import basicConfig, INFO
from os.path import join
from tempfile import TemporaryDirectory
from unittest import main, TestCase

from tips import TotalPartitionFunction


molecule = "H2O"
isotopologue = 1
temperature = 270.
result = 152.19
tolerance = 5
name = "test"


class TestTips(TestCase):

    def test_from_file(self):
        t = TotalPartitionFunction()
        t.download_from_web(molecule)
        with TemporaryDirectory() as directory:
            path = join(directory, "{}.nc".format(name))
            t.write_to_netcdf(path)
            t.read_from_dataset(path)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)

    def test_from_web(self):
        t = TotalPartitionFunction()
        t.download_from_web(molecule)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)

    def test_from_database(self):
        t = TotalPartitionFunction()
        t.download_from_web(molecule)
        with TemporaryDirectory() as directory:
            database = join(directory, "{}.db".format(name))
            t.create_database(database)
            t.load_from_database(molecule, database)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)


if __name__ == "__main__":
    basicConfig(level=INFO)
    main()
