from sqlite3 import connect

from numpy import asarray

from .utils import cm_to_m, download_gfdl_data, GriddedField


class OzoneContinuum(object):
    """Ozone continuum.

    Attributes:
        cross_section: Cross section [cm2].
    """
    def __init__(self, database=None):
        """Initializes object.

        Args:
            database: Path to SQLite database to read from.
        """
        if database is None:
            self.download()
        else:
            self.load_from_database(database)

    def absorption_coefficient(self, grid):
        """Calculates the absorption coefficient.

        Args:
            grid: Numpy array of wavenumbers [cm-1].

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        return self.cross_section.regrid(grid)[:]*cm_to_m*cm_to_m

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            table = "O3_continuum"
            cursor.execute("CREATE TABLE {}(wavenumber REAL, cross_section REAL)".format(table))
            for i in range(self.cross_section.grid.size):
                cursor.execute("""INSERT INTO {}(wavenumber, cross_section)
                                  VALUES (?, ?)""".format(table),
                               (self.cross_section.grid[i], self.cross_section.data[i]))
            connection.commit()

    def download(self):
        """Downloads cross sections from GFDL's FTP site."""
        download_gfdl_data(self, "ozone_continuum", ["ozone_continuum.csv",],
                           ["cross_section",])

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            table = "O3_continuum"
            cursor.execute("SELECT wavenumber, cross_section FROM {}".format(table))
            cross_section = []; grid = []
            for record in cursor.fetchall():
                grid.append(record[0])
                cross_section.append(record[1])
            self.cross_section = GriddedField(asarray(cross_section), asarray(grid))
