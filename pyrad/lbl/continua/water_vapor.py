from sqlite3 import connect

from numpy import asarray, exp

from .utils import cm_to_m, download_gfdl_data, GriddedField, Pa_to_atm
from ...utils.database_utilities import scrub


class _WaterVaporContinuum(object):
    """Water vapor self continuum.

    Attributes:
        c: Broadening coefficient [cm2 atm-1].
        continuum_type: String describing continuum.
        t: Broadening temperature coefficient [K-1].
    """
    def __init__(self, continuum_type, database=None):
        """Initializes object.

        Args:
            continuum_type: String describing continuum ("self" or "foreign").
            database: Path to SQLite database to read from.
        """
        self.continuum_type = continuum_type
        if database is None:
            self.download()
        else:
            self.load_from_database(database)

    def absorption_coefficient(self, temperature, pressure, volume_mixing_ratio, grid):
        """Calculates the absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].
            grid: Numpy array of wavenumbers [cm-1].

        Returns:
            Numpy array of absorption coefficients [m2].
        """
        partial_pressure = self.partial_pressure(pressure, volume_mixing_ratio)
        return (296./temperature)*self.c.regrid(grid)[:]*partial_pressure * \
               exp(self.t.regrid(grid)[:]*(296. - temperature))*cm_to_m*cm_to_m

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            table = "H2O_{}_continuum".format(scrub(self.continuum_type))
            cursor.execute("CREATE TABLE {}(wavenumber REAL, c REAL, t REAL)".format(table))
            for i in range(self.t.grid.size):
                cursor.execute("""INSERT INTO {}(wavenumber, c, t)
                                  VALUES (?, ?, ?)""".format(table),
                               (self.c.grid[i], self.c.data[i], self.t.data[i]))
            connection.commit()

    def download(self, filenames):
        """Downloads water vapor continuum coefficients from GFDL's FTP site.

        Args:
            filenames: List of files to read.
        """
        download_gfdl_data(self, "water_vapor_continuum", filenames, ["c", "t"])

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            table = "H2O_{}_continuum".format(scrub(self.continuum_type))
            cursor.execute("SELECT wavenumber, c, t FROM {}".format(table))
            c = []; grid = []; t = []
            for record in cursor.fetchall():
                grid.append(record[0])
                c.append(record[1])
                t.append(record[2])
            for name, data in zip(["c", "t"], [c, t]):
                setattr(self, name, GriddedField(asarray(data), asarray(grid)))

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        raise NotImplementedError("you must override this function.")


class WaterVaporSelfContinuum(_WaterVaporContinuum):
    """Water vapor self continuum."""
    def __init__(self, database=None):
        """Initializes object.

        Args:
            database: Path to SQLite database to read from.
        """
        super().__init__("self", database)

    def download(self):
        """Downloads water vapor self-continuum coefficients from GFDL's FTP site."""
        super().download(["296MTCKD25_S.csv", "CKDS.csv"])

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        return pressure*volume_mixing_ratio*Pa_to_atm


class WaterVaporForeignContinuum(_WaterVaporContinuum):
    """Water vapor foreign continuum."""
    def __init__(self, database=None):
        """Initializes object.

        Args:
            database: Path to SQLite database to read from.
        """
        super().__init__("foreign", database)

    def download(self):
        """Downloads water vapor continuum foreign-coefficients from GFDL's FTP site."""
        super().download(["296MTCKD25_F.csv", "CKDF.csv"])

    @staticmethod
    def partial_pressure(pressure, volume_mixing_ratio):
        """Calculates the partial pressure.

        Args:
            pressure: Pressure [Pa].
            volume_mixing_ratio: Volume mixing ratio [mol mol-1].

        Returns:
            Partial pressure [atm].
        """
        return pressure*(1. - volume_mixing_ratio)*Pa_to_atm
