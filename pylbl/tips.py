from logging import info
from re import match, search
from sqlite3 import connect
from urllib.request import urlopen

from netCDF4 import Dataset
from numpy import asarray, copy, float32, float64, transpose, searchsorted

from .database_utilities import scrub


tips_reference_temperature = 296.


class MoleculeNotFound(Exception):
    pass


class TotalPartitionFunction(object):
    """Total partition function, using data from TIPS 2017 (doi: 10.1016/j.jqsrt.2017.03.045).

    Attributes:
        data: Numpy array of total partition function values (isotopologue, temperature).
        molecule: Molecule chemical formula.
        temperature: Numpy array of temperatures [K].
    """

    def __init__(self):
        self.data = None
        self.molecule = None
        self.temperature = None

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        connection = connect(database)
        cursor = connection.cursor()
        name = scrub(self.molecule)
        data = transpose(self.data)
        columns = ", ".join(["temperature REAL"] +
                            ["Q_{} REAL".format(i+1) for i in range(data.shape[1])])
        cursor.execute("CREATE TABLE {} ({})".format(name, columns))
        value_subst = ", ".join(["?" for _ in range(data.shape[1] + 1)])
        for i, t in enumerate(self.temperature):
            values = tuple(float64(x) for x in ([t] + data[i, :].tolist()))
            cursor.execute("INSERT INTO {} VALUES ({})".format(name, value_subst), values)
        connection.commit()
        connection.close()

    def download_from_web(self, molecule):
        """Downloads the data from the internet.

        Args:
            molecule: Molecule id.
        """
        url = "http://faculty.uml.edu/Robert_Gamache/Software/temp/Supplementary_file.txt"
        info("Downloading TIPS 2017 data for {} from {}.".format(molecule, url))
        request = urlopen(url)
        records = self.records(request, molecule)
        self.parse_records(records)
        self.molecule = molecule

    @property
    def isotopologue(self):
        return [x for x in range(self.data.shape[0])]

    def load_from_database(self, molecule, database):
        """Loads data from a previously created SQLite database.

        Args:
            molecule: Molecule chemical formula.
            database: Path to SQLite database.
        """
        connection = connect(database)
        cursor = connection.cursor()
        name = scrub(molecule)
        cursor.execute("SELECT * from {}".format(name))
        self.parse_records(cursor.fetchall())
        connection.close()

    def parse_records(self, records):
        """Parses all database records and stores the data.

        Args:
            records: A list of iterable database record values.
        """
        temperature, q = [], []
        for record in records:
            temperature.append(float32(record[0]))
            q.append([float32(x) for x in record[1:]])
        self.temperature = asarray(temperature, dtype=float32)
        self.data = transpose(asarray(q, dtype=float32))
        info("Found data for {} isotopologues at {} temperatures.".format(*self.data.shape))

    def read_from_dataset(self, path):
        """Reads in total partition function table.

        Args:
            path: Path to input file.
        """
        info("Reading TIPS 2017 data from dataset {}.".format(path))
        with Dataset(path, "r") as dataset:
            self.temperature = copy(dataset.variables["temperature"])
            self.data = copy(dataset.variables["total_partition_function"])

    @staticmethod
    def records(request, molecule):
        """Parses the HTTP table for all records related to the input molecule.

        Args:
            request: A urllib.request.Request object.
            molecule: Molecule id.

        Yields:
            A list of floats from a record from the http table.

        Raises:
            MoleculeNotFound: Failed to find the input molecule.
        """
        preamble_byte_length = 42
        request.read(preamble_byte_length)
        record = None
        name_byte_width = 7
        while record != "":
            record = request.read(name_byte_width).decode("utf-8")
            if search(molecule, record):
                eol = ["\0"]
                while eol[-1] != "\n":
                    eol.append(request.read(1).decode("utf-8"))
                num_isotopologues = sum(c == "Q" for c in eol)
                break
            else:
                c = request.read(1).decode("utf-8")
                while c != "\n":
                    c = request.read(1).decode("utf-8")
        else:
            raise MoleculeNotFound("molecule {} not found in TIPS 2017 tables.".format(molecule))
        record_byte_length = 227
        while True:
            record = request.read(record_byte_length).decode("utf-8").strip("\n")
            if not match(r"\s+[0-9]", record):
                return
            yield [float32(x.strip()) for x in record.split()[:(num_isotopologues+1)]]

    def total_partition_function(self, temperature, isotopologue):
        """Interpolates the total partition function values from the TIPS 2017 table.

        Args:
            temperature: Temperature [K].
            isotopologue: Isotopologue id.

        Returns:
            Total partition function.
        """
        i = isotopologue - 1
        j = searchsorted(self.temperature, temperature, side="left") - 1
        return self.data[i,j] + (self.data[i,j+1] - self.data[i,j])*(temperature -
               self.temperature[j])/(self.temperature[j+1] - self.temperature[j])

    def write_to_netcdf(self, path):
        """Writes data to a netCDF dataset.

        Args:
            path: Name of the netCDF4 file that will be created.
        """
        info("Writing TIPS 2017 data to dataset {}.".format(path))
        with Dataset(path, "w") as dataset:
            dataset.createDimension("temperature", self.temperature.size)
            dataset.createDimension("isotopologue", self.data.shape[0])
            v = dataset.createVariable("temperature", float32, dimensions=("temperature",))
            v.setncattr("units", "K")
            v[:] = self.temperature[:]
            v = dataset.createVariable("total_partition_function", float32,
                                       dimensions=("isotopologue", "temperature"))
            v[:,:] = self.data[:,:]