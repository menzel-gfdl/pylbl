from html.parser import HTMLParser
from urllib.request import urlopen


cm_to_m = 0.01 # [m cm-1].


def cross_section_bands(sorted_tables):
    """Organizes a list of sorted Table objects by spectral band.

    Args:
        sorted_tables: A sorted list of Table objects.

    Returns:
        A list of lists of Table objects, organized by spectral band.
    """
    bands = []
    band = [sorted_tables[0],]
    for table in sorted_tables[1:]:
        if table.band_params == band[-1].band_params:
            band.append(table)
        else:
            # Next band.
            bands.append(band)
            band = [table,]
    else:
        bands.append(band)
    return bands


def cross_section_data_files(url):
    """Parses a HITRAN page, looking for absorption cross section data file names.

    Args:
        url: A url string.

    Returns:
        A list of absorption cross section data file names.
    """
    parser = _Parser()
    parser.feed(urlopen(url).read().decode("utf-8"))
    return parser.files


def cross_section_inquiry(tables, id):
    """Handles the case where no HITRAN cross-section files are found.

    Args:
        tables: A list of Table objects.
        id: A string to improve error message clarity.

    Raises:
        A ValueError if no HITRAN absorption cross-section data was located.
    """
    if not tables:
        raise ValueError(" ".join(["Cannot find HITRAN absorption cross-section",
                                   "data files for {}.".format(id)]))


class _Parser(HTMLParser):
    """Parses the collision-induced absorption data file names from the HITRAN website.

    Attributes:
        column: Index of the current column in the HTML table.
        files: List of file names.
        inside_link: Flag telling if a link is currently being processed.
        inside_row: Flag telling if the parser is currently inside a table row.
        inside_table: Flag telling if the parser is currently inside a table.
    """
    def __init__(self):
        HTMLParser.__init__(self)
        self.column = 0
        self.files = []
        self.inside_link = False
        self.inside_row = False
        self.inside_table = False

    def handle_data(self, data):
        """Stores data based on the parser's current state.

        Args:
            data: The data contained in the current HTML tag.
        """
        if self.inside_link:
            self.files.append(data.strip())

    def handle_endtag(self, tag):
        """Resets flags and stores data based on the parser's current state.

        Args:
            tag: Current HTML tag.
        """
        if tag == "table" and self.inside_table:
            self.inside_table = False
        elif tag == "tr" and self.inside_row:
            self.inside_row = False
        elif tag == "a" and self.inside_link:
            self.inside_link = False

    def handle_starttag(self, tag, attrs):
        """Sets flags based on the parser's current state.

        Args:
            tag: Current HTML tag.
            attrs: List of tuples containing HTML tag attributes.
        """
        if tag == "table":
            self.inside_table = True
        elif tag == "tr":
            self.column = 0
            self.inside_row = True
        elif tag == "td" and self.inside_row:
            self.column += 1
        elif tag == "a" and self.inside_row and self.column == 2:
            self.inside_link = True


class Table(object):
    """Container for HITRAN absorption coefficients tables.

    Attributes:
        band_params: Tuple describing spectral band
                     (starting wavenumber [cm-1], ending wavenumber[cm-1], number of points).
        cross_section: Cross sections [cm2 or cm4].
        pressure: Pressure [Pa].
        temperature: Temperature [K].
        wavenumber: Wavenumber [cm-1].
    """
    def __init__(self, temperature, band_params, pressure=None):
        """Creates a Table object.

        Args:
            temperature: Temperature [K].
            band_params: Tuple describing spectral band
                         (starting wavenumber [cm-1], ending wavenumber[cm-1], number of points).
            pressure: Pressure [Torr].
        """
        self.temperature = temperature
        self.band_params = band_params
        if pressure is not None:
            torr_to_Pa = 133.332
            self.pressure = pressure*torr_to_Pa
        self.wavenumber = []
        self.cross_section = []

    def insert(self, wavenumber, cross_section):
        """Appends values to internal lists.

        Args:
            wavenumber: Wavenumber [cm-1].
            cross_section: Cross section [cm2 or cm4].
        """
        self.wavenumber.append(wavenumber)
        self.cross_section.append(cross_section)
