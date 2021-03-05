from logging import getLogger
from re import escape, match, sub
from sqlite3 import connect, IntegrityError
from urllib.request import urlopen

from numpy import abs, argwhere, asarray, linspace, min, searchsorted, zeros
from scipy.interpolate import interp1d, interp2d, LinearNDInterpolator, NearestNDInterpolator

from .utils import cm_to_m, cross_section_bands, cross_section_data_files, \
                   cross_section_inquiry, Table
from ...utils.database_utilities import ascii_table_records, scrub


info = getLogger(__name__).info


class HitranCrossSection(object):
    """HITRAN absorption cross-section database front-end.

    Attributes:
       bands: A list of lists of Table objects containing cross sections.
       molecule: Molecule chemical forumla.
    """
    def __init__(self, molecule, database=None):
        """Downloads HITRAN absorption cross section data tables from the web.

        Args:
            molecule: Molecule chemical formula.
            database: Path to database.
        """
        self.molecule = molecule
        if database is None:
            self.download_from_web()
        else:
            self.load_from_database(database)

    def absorption_coefficient(self, temperature, pressure, grid):
        """Calculates collision-induced absorption coefficients.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [Pa].
            grid: Wavenumber grid [cm-1].

        Returns:
            Array of absorption coefficients [m2].
        """
        cross_section = zeros(grid.size)
        band_left_index = []; band_right_index = []
        for band in self.bands:
            if band[0].band_params[0] > grid[-1] or band[0].band_params[1] < grid[0]:
                # This band is outside the input grid, so skip it.
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
                # Only one temperature and pressure, so use 1D interpolation.
                x = asarray(band[0].wavenumber)
                y = asarray(band[0].cross_section)
                cross_section[left:right] = interp1d(x, y)(grid[left:right])
            else:
                temperatures = [x.temperature for x in band]
                pressures = [x.pressure for x in band]
                x = asarray(band[0].wavenumber)
                z = asarray([x.cross_section for x in band])
                if len(set(temperatures)) == 1:
                    # Only one temperature, so use 2D interpolation.
                    y = asarray(pressures)
                    cross_section[left:right] = interp2d(x, y, z)(grid[left:right], pressure)
                elif len(set(pressures)) == 1:
                    # Only one pressure, so use 2D interpolation.
                    y = asarray(temperatures)
                    cross_section[left:right] = interp2d(x, y, z)(grid[left:right], temperature)
                else:
                    # Multiple temperatures/pressures, so use ND interpolation.  For each
                    # wavenumber triangulate in temperature and pressure space.
                    x = asarray(temperatures)
                    y = asarray(pressures)
                    data = zeros(z.shape[-1])
                    for i in range(z.shape[-1]):
                        if z.shape[0] >= 4:
                            f = LinearNDInterpolator(list(zip(x, y)), z[:,i], fill_value=0.)
                        else:
                            # If there isn't enough data (4 temperature, pressure points),
                            # then default to nearest-neighbor interpolation.
                            f = NearestNDInterpolator(list(zip(x, y)), z[:,i])
                        data[i] = f(temperature, pressure)
                    # Interpolate to the input spectral grid.
                    cross_section[left:right] = interp1d(asarray(band[0].wavenumber),
                                                         data)(grid[left:right])
        # Remove non-physical negative values.
        cross_section[argwhere(cross_section < 0.)] = 0.
        return cross_section*cm_to_m*cm_to_m

    def download_from_web(self):
        """Downloads HITRAN collision-induced absorption cross-sections from the web."""
        url = "http://hitran.org/data/xsec"
        tables = []; temperatures = []; pressures = []; band_params = []
        for f in cross_section_data_files(url):
            if match(escape(self.molecule) + r"_.*" + escape(".xsc"), f):
                address = "{}/{}".format(url, f)
                info("Downloading absorption cross sections from {}.".format(address))
                table = self.parse_records(self.records(urlopen(address)))
                if tables and table.temperature in temperatures and table.band_params in \
                   band_params and table.pressure in pressures:
                    # Skip duplicate data to prevent nans when interpolating.
                    continue
                temperatures.append(table.temperature)
                band_params.append(table.band_params)
                pressures.append(table.pressure)
                tables.append(table)
        self.bands = self.create_bands(tables, self.molecule)

    @staticmethod
    def create_bands(tables, molecule):
        """Organize the cross-section tables into bands.

        Args:
            tables: List of Table objects.
            molecule: String chemical formula.

        Returns:
            A list of lists of Table objects containing cross sections.
        """
        cross_section_inquiry(tables, molecule)
        sorted_tables= sorted(sorted(sorted(sorted(tables, key=lambda x: x.pressure),
                                            key=lambda x: x.temperature),
                                      key=lambda x: x.wavenumber[-1]),
                              key=lambda x: x.wavenumber[0])
        return cross_section_bands(sorted_tables)

    def parse_records(self, records):
        """Parses all database records and stores the data.

        Args:
            records: A list of list containing database record values.

        Returns:
            A Table object.
        """
        torr_to_Pa = 133.332
        for i, record in enumerate(records):
            if i == 0:
                temperature, pressure = (float(x.strip()) for x in record[4:6])
                w0, wn, n = (float(x.strip()) for x in record[1:4])
                band_params = tuple([w0, wn, int(n)])
                table = Table(temperature, band_params, pressure*torr_to_Pa)
                table.wavenumber = linspace(w0, wn, int(n), endpoint=True)
            else:
                table.cross_section += [float(x) for x in record]
        if len(table.cross_section) > len(table.wavenumber):
            # For some reason I cannot understand, HITRAN sometimes adds extra zeros
            # at the end of a file.  Ignore these.
            table.cross_section = table.cross_section[:len(table.wavenumber)]
        return table

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
        molecule = sub("-", "_", scrub(self.molecule))
        with connect(database) as connection:
            cursor = connection.cursor()
            cursor.execute("PRAGMA foreign_keys = ON")

            # Create a table to hold band parameters.
            band_table = "{}_cross_section_band_parameters".format(molecule)
            cursor.execute("""CREATE TABLE {}(band_id INTEGER PRIMARY KEY, lower_bound REAL,
                              upper_bound REAL, size INTEGER)""".format(band_table))
            for band in self.bands:
                cursor.execute("""INSERT INTO {}(lower_bound, upper_bound, size)
                                  VALUES (?, ?, ?)""".format(band_table),
                               band[0].band_params)
            connection.commit()

            # Create a table to hold temperatures.
            temperature_table = "{}_cross_section_temperatures".format(molecule)
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

            # Create a table to hold pressures.
            pressure_table = "{}_cross_section_pressures".format(molecule)
            cursor.execute("""CREATE TABLE {}(pressure_id INTEGER PRIMARY KEY,
                              pressure REAL UNIQUE)""".format(pressure_table))
            for band in self.bands:
                for table in band:
                    try:
                        cursor.execute("""INSERT INTO {}(pressure)
                                          VALUES (?)""".format(pressure_table),
                                       (table.pressure,))
                    except IntegrityError:
                        # Prevent duplicate pressure entries.
                        continue
            connection.commit()

            # Create a table of cross sections.
            cross_section_table = "{}_cross_sections".format(molecule)
            cursor.execute("""CREATE TABLE {}(wavenumber REAL, cross_section REAL,
                              band INTEGER REFERENCES {}(band_id), temperature INTEGER
                              REFERENCES {}(temperature_id), pressure INTEGER REFERENCES
                              {}(pressure_id))""".format(cross_section_table, band_table,
                                                         temperature_table, pressure_table))
            for band in self.bands:
                cursor.execute("""SELECT band_id FROM {} WHERE lower_bound == {} AND
                                  upper_bound == {} AND size == {}""".format(band_table,
                                                                             *band[0].band_params))
                band_id = cursor.fetchone()[0]
                for table in band:
                    cursor.execute("""SELECT temperature_id FROM {} WHERE temperature ==
                                      {}""".format(temperature_table, table.temperature))
                    temperature_id = cursor.fetchone()[0]
                    cursor.execute("""SELECT pressure_id FROM {} WHERE pressure ==
                                      {}""".format(pressure_table, table.pressure))
                    pressure_id = cursor.fetchone()[0]
                    for wavenumber, cross_section in zip(table.wavenumber, table.cross_section):
                        cursor.execute("""INSERT INTO {} VALUES
                                          (?, ?, ?, ?, ?)""".format(cross_section_table),
                                       (wavenumber, cross_section, band_id, temperature_id,
                                        pressure_id))
            connection.commit()

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database.
        """
        molecule = sub("-", "_", scrub(self.molecule))
        with connect(database) as connection:
            cursor = connection.cursor()
            cursor.execute("PRAGMA foreign_keys = ON")

            # Read the band parameters.
            name = "{}_cross_section_band_parameters".format(molecule)
            cursor.execute("SELECT band_id, lower_bound, upper_bound, size FROM {}".format(name))
            band_params = [record for record in cursor.fetchall()]

            # Read records from the database table.
            name = "{}_cross_sections".format(molecule)
            tables = []
            for band in band_params:
                cursor.execute("""SELECT wavenumber, cross_section, temperature, pressure
                                  FROM {} WHERE band == {}""".format(name, band[0]))
                temperature = None; pressure = None
                for record in cursor.fetchall():
                    if record[2] == temperature and record[3] == pressure:
                        tables[-1].insert(*record[:2])
                    else:
                        # Get the temperature and pressure.
                        foreign_name = "{}_cross_section_temperatures".format(molecule)
                        cursor.execute("""SELECT temperature FROM {} WHERE temperature_id
                                          == {}""".format(foreign_name, record[2]))
                        t = cursor.fetchone()[0]
                        foreign_name = "{}_cross_section_pressures".format(molecule)
                        cursor.execute("""SELECT pressure FROM {} WHERE pressure_id
                                          == {}""".format(foreign_name, record[3]))
                        p = cursor.fetchone()[0]
                        tables.append(Table(t, band[1:], p))
                        tables[-1].insert(*record[:2])
                        temperature, pressure = record[2:]
        self.bands = self.create_bands(tables, molecule)
