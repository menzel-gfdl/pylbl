from logging import basicConfig, INFO
from tempfile import NamedTemporaryFile
from unittest import main, TestCase

from pyrad.lbl.tips import TotalPartitionFunction


molecule = "H2O"
isotopologue = 1
temperature = 270.
result = 152.19
tolerance = 5


class TestTips(TestCase):

    def test_from_file(self):
        t = TotalPartitionFunction(molecule)
        with NamedTemporaryFile(suffix=".nc") as path:
            t.write_to_netcdf(path.name)
            t.read_from_dataset(path.name)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)

    def test_from_web(self):
        t = TotalPartitionFunction(molecule)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)

    def test_from_database(self):
        t = TotalPartitionFunction(molecule)
        with NamedTemporaryFile(suffix=".db") as database:
            t.create_database(database.name)
            t.load_from_database(database.name)
        q = t.total_partition_function(temperature, isotopologue)
        self.assertAlmostEqual(q, result, places=tolerance)


if __name__ == "__main__":
    basicConfig(level=INFO)
    main()
