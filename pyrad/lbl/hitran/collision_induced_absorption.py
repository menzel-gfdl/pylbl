from logging import getLogger
from re import escape, match
from sqlite3 import connect, IntegrityError
from urllib.request import urlopen

from numpy import abs, argwhere, asarray, float64, min, searchsorted, zeros
from scipy.interpolate import interp1d, interp2d

from .utils import cm_to_m, cross_section_bands, cross_section_data_files, \
                   cross_section_inquiry, Table
from ...utils.database_utilities import ascii_table_records, scrub


info = getLogger(__name__).info


class HitranCIA(object):
    """HITRAN collision-inducd absorption database front-end.

    Attributes:
        bands: A list of lists of Table objects containing cross sections.
        broadener: Molecule chemical formula for the broadener.
        molecule: Molecule chemical forumla.
    """

    def __init__(self, molecule, broadener, database=None):
        """Downloads HITRAN collision-induced absorption cross sections from the web.

        Args:
            molecule: Molecule chemical forumla.
            broadener: Molecule chemical formula for the broadener.
            database: Path to database.
        """
        self.broadener = broadener
        self.molecule = molecule
        if database is None:
            self.download_from_web()
        else:
            self.load_from_database(database)

    def absorption_coefficient(self, temperature, grid):
        """Calculates collision-induced absorption coefficients.

        Args:
            temperature: Temperature [K].
            grid: Spectral grid array [cm-1].

        Returns:
            Array of absorption coefficients [m5].
        """
        cross_section = zeros(grid.size)
        band_left_index, band_right_index = [], []
        for band in self.bands:
            if band[0].band_params[0] > grid[-1] or band[0].band_params[1] < grid[0]:
                # This band is outside of the input grid, so skip it.
                continue

            # Find part of the input grid that overlaps with this spectral band.
            left = searchsorted(grid, band[0].band_params[0])
            right = searchsorted(grid, band[0].band_params[1])

            # Check for spectral overlap between bands.
            override = True
            for i, (lower, upper) in enumerate(zip(band_left_index, band_right_index)):
                if left <= upper:
                    # There is some band overlap.
                    if right <= upper:
                        # This band lies entirely inside another band.
                        t_self = min(abs(asarray([x.temperature for x in band]) -
                                         temperature))
                        t_other = min(abs(asarray([x.temperature for x in self.bands[i]]) -
                                          temperature))
                        override = t_self < t_other
                    else:
                        # This band extends beyond the band it overlaps with.
                        left = upper
                    break
            if not override:
                # Skip this band because the broader band's cross sections are at
                # a closer temperature.
                continue
            band_left_index.append(left)
            band_right_index.append(right)
            if len(band) == 1:
                # Only one temperature, so use 1D interpolation.
                x = asarray(band[0].wavenumber)
                y = asarray(band[0].cross_section)
                cross_section[left:right] = interp1d(x, y)(grid[left:right])
            else:
                # Multiple temperatures, so use 2D interpolation.
                x = asarray(band[0].wavenumber)
                y = asarray([x.temperature for x in band])
                z = asarray([x.cross_section for x in band])
                cross_section[left:right] = interp2d(x, y, z)(grid[left:right], temperature)
        # Remove non-physical negative values.
        cross_section[argwhere(cross_section < 0.)] = 0.
        return cross_section*cm_to_m*cm_to_m*cm_to_m*cm_to_m*cm_to_m

    def download_from_web(self):
        """Downloads HITRAN collision-induced absorption cross-sections from the web."""
        urls = ["http://hitran.org/data/CIA",]
        interaction = "{}-{}".format(self.molecule, self.broadener)
        tables = []; temperatures = []; band_params = []
        for url in urls:
            for f in cross_section_data_files(url):
                if match(escape(interaction) + r"[A-Za-z0-9_]*" + escape(".cia"), f):
                    address = "{}/{}".format(url, f)
                    info("Downloading collision-absorption cross sections from {}.".format(address))
                    for table in self.parse_records(self.records(urlopen(address))):
                        if tables and table.temperature in temperatures and \
                           table.band_params in band_params:
                            # Skip duplicate data to prevent nans when interpolating.
                            continue
                        temperatures.append(table.temperature)
                        band_params.append(table.band_params)
                        tables.append(table)
        self.bands = self.create_bands(tables, interaction)

    @staticmethod
    def create_bands(tables, interaction):
        """Organize the cross-section tables into bands.

        Args:
            tables: List of Table objects.
            interaction: String describing the collision.

        Returns:
            A list of lists of Table objects containing cross sections.
        """
        cross_section_inquiry(tables, interaction)
        sorted_tables = sorted(sorted(sorted(tables, key=lambda x: x.temperature),
                                      key=lambda x: x.wavenumber[-1]),
                               key=lambda x: x.wavenumber[0])
        return cross_section_bands(sorted_tables)

    def parse_records(self, records):
        """Parses all database records and stores the data.

        Args:
            records: A list of list containing database record values.

        Returns:
            A list of Table objects.

        Raises:
            A ValueError if a record has an unexpected format.
        """
        tables = []
        for record in records:
            if len(record) > 4:
                temperature = float(record[4].strip())
                band_params = tuple([float(x.strip()) for x in record[1:3]] +
                                    [int(record[3].strip()),])
                tables.append(Table(temperature, band_params))
            elif len(record) == 2 or len(record) == 3:
                tables[-1].insert(*[float(x.strip()) for x in record[:2]])
            else:
                raise ValueError("record does not meet expected format:\n{}".format(record))
        return tables

    def records(self, response):
        """Parses the HTTP table for all records related to the input molecule.

        Args:
            response: A http.client.HTTPResponse object.

        Yields:
            A list values from a record from the http table.
        """
        for line in ascii_table_records(response):
            yield line.split()

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        interaction = scrub("{}_{}".format(self.molecule, self.broadener))
        with connect(database) as connection:
            cursor = connection.cursor()
            cursor.execute("PRAGMA foreign_keys = ON")

            # Create a table to hold band parameters.
            band_table = "{}_cia_band_parameters".format(interaction)
            cursor.execute("""CREATE TABLE {}(band_id INTEGER PRIMARY KEY, lower_bound REAL,
                              upper_bound REAL, size INTEGER)""".format(band_table))
            for band in self.bands:
                cursor.execute("""INSERT INTO {}(lower_bound, upper_bound, size)
                                  VALUES (?, ?, ?)""".format(band_table),
                               band[0].band_params)
            connection.commit()

            # Create a table to hold temperatures.
            temperature_table = "{}_cia_temperatures".format(interaction)
            cursor.execute("""CREATE TABLE {}(temperature_id INTEGER PRIMARY KEY,
                              temperature REAL UNIQUE)""".format(temperature_table))
            for band in self.bands:
                for table in band:
                    try:
                        cursor.execute("""INSERT INTO {}(temperature)
                                          VALUES (?)""".format(temperature_table),
                                       (table.temperature,))
                    except IntegrityError:
                        # Prevent duplicate temperature entries.
                        continue
            connection.commit()

            # Create a table of cross sections.
            cross_section_table = "{}_cia_cross_sections".format(interaction)
            cursor.execute("""CREATE TABLE {}(wavenumber REAL, cross_section REAL,
                              band INTEGER REFERENCES {}(band_id), temperature INTEGER
                              REFERENCES {}(temperature_id))""".format(cross_section_table,
                                                                       band_table,
                                                                       temperature_table))
            for band in self.bands:
                cursor.execute("""SELECT band_id FROM {} WHERE lower_bound == {} AND
                                  upper_bound == {} AND size == {}""".format(band_table,
                                                                             *band[0].band_params))
                band_id = cursor.fetchone()[0]
                for table in band:
                    cursor.execute("""SELECT temperature_id FROM {} WHERE temperature ==
                                      {}""".format(temperature_table, table.temperature))
                    temperature_id = cursor.fetchone()[0]
                    for wavenumber, cross_section in zip(table.wavenumber, table.cross_section):
                        cursor.execute("""INSERT INTO {} VALUES
                                          (?, ?, ?, ?)""".format(cross_section_table),
                                       (wavenumber, cross_section, band_id, temperature_id))
            connection.commit()

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database.
        """
        interaction = scrub("{}_{}".format(self.molecule, self.broadener))
        with connect(database) as connection:
            cursor = connection.cursor()
            cursor.execute("PRAGMA foreign_keys = ON")

            # Read the band parameters.
            name = "{}_cia_band_parameters".format(interaction)
            cursor.execute("SELECT band_id, lower_bound, upper_bound, size FROM {}".format(name))
            band_params = [record for record in cursor.fetchall()]

            # Read the temperatures.
            name = "{}_cia_temperatures".format(interaction)
            cursor.execute("SELECT temperature_id, temperature FROM {}".format(name))
            temperatures = [record for record in cursor.fetchall()]

            # Read records from the database table.
            name = "{}_cia_cross_sections".format(interaction)
            tables = []
            for band in band_params:
                for temperature in temperatures:
                    cursor.execute("""SELECT wavenumber, cross_section FROM {} WHERE band ==
                                    {} AND temperature == {}""".format(name, band[0],
                                                                          temperature[0]))
                    wavenumber = []; cross_section = []
                    for record in cursor.fetchall():
                        wavenumber.append(record[0])
                        cross_section.append(record[1])
                    if wavenumber:
                        tables.append(Table(temperature[1], band[1:]))
                        for w, xsec in zip(wavenumber, cross_section):
                            tables[-1].insert(w, xsec)
        self.bands = self.create_bands(tables, interaction)
