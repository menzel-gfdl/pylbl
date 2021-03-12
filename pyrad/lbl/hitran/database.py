from collections import OrderedDict
from logging import getLogger
from sqlite3 import connect
from urllib.request import urlopen

from numpy import asarray

from .isotopologues import isotopologues
from .line_parameters import linear_molecule_quantum_numbers, linear_molecule_state_vars, \
                             PARAMETERS
from .molecules import molecules
from .spectral_lines import SpectralLines
from ...utils.database_utilities import ascii_table_records, scrub, SQL_TYPES


info = getLogger(__name__).info
warning = getLogger(__name__).warning


class Hitran(object):
    """HITRAN database front-end.

    Attributes:
       isotopologues: Lists of Isotopologue objects.
       line_profile: Doppler, Lorentz, or Voigt object.
       molecule: Molecule chemical forumla.
       molecule_id: HITRAN molecule id.
       parameters: List of HitranParameter objects.
    """

    def __init__(self, molecule, line_profile, isotopologue=None, database=None):
        """Creates dictionaries of HITRAN ids by scraping the HITRAN website."""
        self.molecule = molecule
        self.line_profile = line_profile
        base_parameters = ["id", "iso", "center", "strength", "elower", "delta_air"]
        all_parameters = base_parameters + \
                         [x for x in line_profile.parameters if x not in base_parameters]
        self.parameters = [PARAMETERS[x] for x in all_parameters]
        info(" ".join(["Using parameters"] + [x.shortname for x in self.parameters]))
        self.molecule_id = molecules("https://hitran.org/docs/molec-meta/")[molecule]
        self.isotopologues = isotopologue if isotopologue is not None else \
            isotopologues("https://hitran.org/docs/iso-meta/")[molecule]
        for parameter in self.parameters:
            setattr(self, parameter.shortname, list())
        if database is None:
            self.download_from_web()
        else:
            self.load_from_database(database)

    def create_database(self, database):
        """Creates/ingests data into a SQLite database.

        Args:
            database: Path to SQLite database that will be create/added to.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            name = "{}_lines".format(scrub(self.molecule))
            shortnames, types = [], []
            for x in self.parameters:
                if x.dtype is linear_molecule_quantum_numbers:
                    q = getattr(self, x.shortname)[0]
                    types += [SQL_TYPES[i] for i in linear_molecule_state_vars.values()]
                    state = "upper" if x.api_name == "statep" else "lower"
                    shortnames += ["{}_{}".format(i, state) for i in
                                   linear_molecule_state_vars.keys()]
                else:
                    types.append(SQL_TYPES[x.dtype])
                    shortnames.append(x.shortname)
            print(shortnames)
            value_subst = ", ".join(["?" for _ in shortnames])
            header = ", ".join(["{} {}".format(shortname, type_str)
                                for shortname, type_str in zip(shortnames, types)])
            cursor.execute("CREATE TABLE {} ({})".format(name, header))
            for i in range(getattr(self, self.parameters[0].shortname).shape[0]):
                # Extra type cast is needed because sqlite3 cannot handle numpy int
                # objects as INTEGER values.
                values = []
                for x in self.parameters:
                    if x.dtype is linear_molecule_quantum_numbers:
                        values += [i(q) for i, q in zip(linear_molecule_state_vars.values(),
                                                        getattr(self, x.shortname)[i].values())]
                    else:
                        try:
                            values.append(x.dtype(getattr(self, x.shortname)[i]))
                        except IndexError:
                            print(x.shortname)
                            raise
                cursor.execute("INSERT INTO {} VALUES ({})".format(name, value_subst),
                               tuple(values))
            connection.commit()

    def download_from_web(self, lower_bound=0., upper_bound=10.e6):
        """Downloads HITRAN molecular line parameters from the web.

        Args:
            lower_bound: Lower bound of spectral range [cm-1], inclusive.
            upper_bound: Upper bound of spectral range [cm-1], inclusive.
        """
        options = [("iso_ids_list", ",".join([str(x.id) for x in self.isotopologues])),
                   ("numin", lower_bound),
                   ("numax", upper_bound),
                   ("request_params", ",".join([x.api_name for x in self.parameters]))]
        html_options = "&".join(["{}={}".format(*x) for x in options])
        url = "http://hitran.org/lbl/api?{}".format(html_options)
        info("Downloading Hitran database for {} from {}.".format(self.molecule, url))
        self.parse_records(self.records(urlopen(url)))

    def spectral_lines(self, total_partition_function):
        return SpectralLines(self, total_partition_function)

    def load_from_database(self, database):
        """Loads data from a previously created SQLite database.

        Args:
            database: Path to SQLite database.
        """
        with connect(database) as connection:
            cursor = connection.cursor()
            name = "{}_lines".format(scrub(self.molecule))
            shortnames = [x.shortname for x in self.parameters
                          if x.shortname != "q1" and x.shortname != "q2"]
            quantum_numbers = len(shortnames) != len(self.parameters)
            if quantum_numbers:
                #Quantum numbers needed.
                all_parameters = self.parameters.copy()
                self.parameters = [x for x in all_parameters if x.shortname in shortnames]
            #Read all the non-quantum number parameters out of the database.
            columns = ", ".join(shortnames)
            cursor.execute("SELECT {} from {}".format(columns, name))
            self.parse_records(cursor.fetchall())
            if quantum_numbers:
                #Now read in the quantum numbers and create the ordered dictionaries.
                for q, state in zip(["q1", "q2"], ["upper", "lower"]):
                    data = []
                    shortnames = [x for x in linear_molecule_state_vars.keys()]
                    types = [x for x in linear_molecule_state_vars.values()]
                    columns = ", ".join(["{}_{}".format(x, state) for x in shortnames])
                    cursor.execute("SELECT {} from {}".format(columns, name))
                    for record in cursor.fetchall():
                        data.append(OrderedDict((x.lower(), y(z)) for x, y, z in
                                                zip(shortnames, types, record)))
                    setattr(self, q, asarray(data))
                self.parameters = all_parameters

    def parse_records(self, records):
        """Parses all database records and stores the data.

        Args:
            records: A list of dictionaries containing database record values.
        """
        for record in records:
            data = []
            for x, y in zip(self.parameters, record):
                try:
                    data.append(x.dtype(y))
                except ValueError as e:
                    if "#" in str(e):
                        if x.shortname == "n_self":
                            data.append(None)
                        else:
                            warning("bad data value in database record:\n{}".format(record))
                            continue
                    else:
                        raise
            if len(data) < len(self.parameters):
                continue
            for x, datum in zip(self.parameters, data):
                self.__dict__[x.shortname].append(datum)
            if "n_self" in self.__dict__ and self.n_self[-1] is None:
                self.n_self[-1] = self.n_air[-1]
        for x in self.parameters:
            setattr(self, x.shortname, asarray(getattr(self, x.shortname)))
        info("Found data for {} lines for {}.".format(getattr(self, self.parameters[0].shortname).shape[0],
                                                      self.molecule))

    def records(self, response):
        """Parses the HTTP table for all records related to the input molecule.

        Args:
            response: A http.client.HTTPResponse object.

        Yields:
            A dictionary of values from a record from the http table.
        """
        for line in ascii_table_records(response):
            yield line.split(",")
