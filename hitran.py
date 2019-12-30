from collections import namedtuple, OrderedDict
from re import match
from sqlite3 import connect
from urllib.request import urlopen

from numpy import argsort, asarray, exp, finfo, float32, float64, log, power, searchsorted, sqrt, sum, zeros

from database_utilities import scrub
from hitran_html_parsers import HitranIsotopologueHTMLParser, HitranMoleculeIdHTMLParser
from line_profiles import Doppler, Lorentz, Voigt
from tips import tips_reference_temperature


RecordFormat = namedtuple("RecordFormat", ["position", "width", "name", "type"])
columns = (RecordFormat(position=0, width=2, name="hitran_id", type=int),
           RecordFormat(position=2, width=1, name="isotopologue", type=int),
           RecordFormat(position=3, width=12, name="center", type=float32),
           RecordFormat(position=15, width=10, name="strength", type=float32),
           RecordFormat(position=35, width=5, name="air_broadened_halfwidth", type=float32),
           RecordFormat(position=40, width=5, name="self_broadened_halfwidth", type=float32),
           RecordFormat(position=45, width=10, name="lower_state_energy", type=float32),
           RecordFormat(position=55, width=4, name="air_broadened_temperature_dependence", type=float32),
           RecordFormat(position=59, width=8, name="air_broadened_pressure_shift", type=float32))
members = OrderedDict([("isotopologue", "iso"), ("center", "v"), ("strength", "s"),
                       ("air_broadened_halfwidth", "gamma_air"),
                       ("self_broadened_halfwidth", "gamma_self"), ("lower_state_energy", "en"),
                       ("air_broadened_temperature_dependence", "n"),
                       ("air_broadened_pressure_shift", "delta")])


class RecordError(BaseException):
    """Error detected while parsing HITRAN spectral line record."""
    pass


class EmptySpectraError(BaseException):
    """No spectral lines found for input molecule/spectral range."""
    pass


class HitranSpectralLine(object):
    """HITRAN spectral parameters for a single molecular line.

    Attributes:
        delta: Air-broadended pressure shift [cm-1 atm-1].
        en: Transition lower state energy [cm-1].
        gamma_air: Air-broadened halfwidth [cm-1 atm-1].
        gamma_self: Self-broadened halfwidth [cm-1 atm-1].
        iso: HITRAN isotopologue id.
        n: Air-broadened temperature dependence power.
        s: Line strength [cm-1].
        v: Transition wavenumber [cm-1].
    """
    def __init__(self):
        for column in columns[1:]:
            setattr(self, members[column.name], None)

    def parse_ascii_record(self, record, id, lower_bound, upper_bound):
        """Parses spectral line parameters from a HITRAN database record.

        Args:
            record: HITRAN database record.
            id: HITRAN molecule id.
            lower_bound: Spectral lower bound [cm-1], inclusive.
            upper_bound: Spectral upper bound [cm-1], inclusive.

        Raises:
            ValueError: An invalid isotopologue id was detected.
            RecordError: A non-matching molecule id or out-of-range transition frequency was detected.
        """
        for column in columns:
            try:
                value = column.type(record[column.position:(column.position + column.width)])
            except ValueError as e:
                letter = match(r"invalid literal for int\(\) with base 10: '([A-Za-z])'", str(e))
                if column.name == "isotopologue" and letter:
                    value = ord(letter.group(1).upper()) - ord("A") + 11
                else:
                    raise
            if column.name == "hitran_id":
                if value == id:
                    continue
                else:
                    raise RecordError("non-matching molecule id.")
            elif column.name == "isotopologue" and value == 0:
                value = 10
            elif column.name == "center" and (value < lower_bound or value > upper_bound):
                raise RecordError("line center not in input spectral range.")
            setattr(self, members[column.name], value)


class HitranSpectralLines(object):
    """Hitran spectral line parameters.

    Attributes:
        delta: Numpy array of air-broadended pressure shifts [cm-1 atm-1] (lines).
        en: Numpy array of transition lower state energies [cm-1] (lines).
        gamma_air: Numpy array of air-broadened halfwidths [cm-1 atm-1] (lines).
        gamma_self: Numpy array of self-broadened halfwidths [cm-1 atm-1] (lines).
        iso: Numpy array of HITRAN isotopologue ids (lines).
        mass: Numpy array of isotopologue masses (lines).
        n: Numpy array of air-broadened temperature dependence powers (lines).
        s: Numpy array of line strengths [cm-1] (lines).
        q: TotalPartitionFunction object.
        t_ref: HITRAN reference temperature.
        v: Numpy array of transition wavenumbers [cm-1] (lines).
    """

    def __init__(self, line_parameters, isotopologues, total_partition_function):
        """Sorts line parameters by transition wavenumber, and partially corrects the line strengths.

        Args:
            line_parameters: List of HitranSpectralLine objects.
            isotopologues: List of isotopologue namedtuples.
            total_partition_function: TotalPartitionFunction object.

        Raises:
            EmptySpectraError: No molecular line parameters are detected.
        """
        if not line_parameters:
            raise EmptySpectraError("no spectral lines in the input spectral range.")
        num_lines = len(line_parameters)
        parameters = sorted(line_parameters, key=lambda line: line.v)
        for key, value in vars(parameters[0]).items():
            setattr(self, key, zeros(num_lines, dtype=type(value)))
        for i, p in enumerate(parameters):
            for key, value in vars(p).items():
                #Convert from array-of-structs to struct-of-arrays.
                self.__dict__[key][i] = value
        self.q = total_partition_function
        self.s[:] *= self.temperature_correct_line_strength(self.q, tips_reference_temperature,
                                                            self.iso[:], self.en[:], self.v[:])
        self.mass = asarray([isotopologues[x-1].mass for x in self.iso])

    def absorption_coefficient(self, temperature, pressure, partial_pressure, wavenumber,
                               line_profile, cut_off=25.):
        """Calculates the absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
            wavenumber: Numpy array of wavenumbers [cm-1] (wavenumber).
            line_profile: Function that calculates the line profile [cm-1].
            cut_off: Distance [cm-1] from the transition frequency where the line is cut off.

        Returns:
            Numpy array of absorption coefficients [cm-2] (wavenumber).
        """
        s = self.correct_line_strengths(temperature)
        vstar = self.pressure_shift_transition_wavenumbers(pressure)

        #Sort lines in order of transition wavenumber.
        inds = vstar.argsort()
        vstar_sorted = vstar[inds]
        s_sorted = s[inds]

        if isinstance(line_profile, Lorentz) or isinstance(line_profile, Voigt):
            gamma = self.pressure_broadened_halfwidths(pressure, partial_pressure, temperature)
            gamma_sorted = gamma[inds]
        if isinstance(line_profile, Doppler) or isinstance(line_profile, Voigt):
            mass_sorted = self.mass[inds]
            alpha = self.doppler_broadened_halfwidth(temperature, mass_sorted, vstar_sorted)

        k = zeros(wavenumber.size)
        for i, v in enumerate(wavenumber):
            l = searchsorted(vstar_sorted, v - cut_off, side="left")
            r = searchsorted(vstar_sorted, v + cut_off, side="right")
            line_profile.gamma = gamma_sorted[l:r]
            line_profile.alpha = alpha[l:r]
            k[i] = sum(s_sorted[l:r]*line_profile.profile(v, vstar_sorted[l:r]))
        return k

    def correct_line_strengths(self, temperature):
        """Temperature-corrects the line strengths.

        Args:
            temperature: Temperature [K].

        Returns:
            Numpy array of corrected line strengths [cm-1] (lines).
        """
        return self.s[:]*(1./self.temperature_correct_line_strength(self.q, temperature,
                                                                    self.iso[:], self.en[:],
                                                                    self.v[:]))

    @staticmethod
    def doppler_broadened_halfwidth(temperature, mass, transition_wavenumber):
        """Calculate the doppler-broadened line halfwidth.

        Args:
            temperature: Temperature [K].
            mass: Molecular mass [g].
            transition_wavenumber: Transition wavenumber [cm-1].

        Returns:
            Doppler-broadened line halfwidth [cm-1].
        """
        m = mass/6.023e23 #Mass/Avagadro's number [g].
        c = 2.99792458e10 #Speed of light [cm s-1].
        kb = 1.380658e-16 #Boltzmann constant [erg K-1].
        return sqrt(log(2))*transition_wavenumber*sqrt(2*kb*temperature/(m*c*c))

    def pressure_broadened_halfwidths(self, pressure, partial_pressure, temperature):
        """Calculates pressure-broadened line halfwidths.

        Args:
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
            temperature: Temperature [K]

        Returns:
            Numpy array of pressure-broadened line halfwidths [cm-1] (lines).
        """
        return power((296./temperature), self.n[:])*(self.gamma_air[:]*
               (pressure - partial_pressure) + self.gamma_self[:]*partial_pressure)

    def pressure_shift_transition_wavenumbers(self, pressure):
        """Pressure-shifts transition wavenumbers.

        Args:
            pressure: Pressure [atm].

        Returns:
            Numpy array of pressure-shifted transition wavenumbers [cm-1] (lines).
        """
        return self.v[:] + self.delta[:]*pressure

    @staticmethod
    def temperature_correct_line_strength(q, t, iso, en, v):
        """Temperature-corrects a line strength.
        Args:
            q: TotalPartitionFunction object.
            t: Temperature [K].
            iso: HITRAN isotopologue id.
            en: Transition lower state energy [cm-1].
            v: Transition wavenumber [cm-1].

        Returns:
            Corrected line strengths [cm-1].
        """
        c2 = -1.4387686
        #Divide-by-zeros may occur for transition wavenumbers close to zero, like those
        #for the O16-O17 isotopologue of O2.
        return q.total_partition_function(t, iso)/(exp(c2*en/t)*(1. - exp(c2*v/t)))


class HitranDatabase(object):
    """HITRAN database front-end.

    Attributes:
       molecule_ids: Dictionary mapping chemical formulae to HITRAN molecule ids.
       isotopologues: Dictionary mapping chemical formuale to lists of Isotopologue objects.
    """

    def __init__(self):
        """Creates dictionaries of HITRAN ids by scraping the HITRAN website."""
        parser = HitranMoleculeIdHTMLParser()
        url = "https://hitran.org/docs/molec-meta/"
        parser.feed(urlopen(url).read().decode("utf-8"))
        self.molecule_ids = parser.molecule_ids
        parser = HitranIsotopologueHTMLParser()
        url = "https://hitran.org/docs/iso-meta/"
        parser.feed(urlopen(url).read().decode("utf-8"))
        self.isotopologues = parser.isotopologues

    def download_from_web(self, molecule, isotopologues="all", lower_bound=0.,
                          upper_bound=finfo("f4").max):
        """Downloads HITRAN molecular line parameters from the web.

        Args:
            molecule: String molecule chemical formula (e.g., "H2O").
            isotopologues: List of HITRAN isotopologue ids.
            lower_bound: Lower bound of spectral range [cm-1], inclusive.
            upper_bound: Upper bound of spectral range [cm-1], inclusive.
        """
        url = "http://hitran.org/lbl/api?"
        if isotopologues == "all":
            iso = ",".join([str(x.id) for x in self.isotopologues[molecule]])
        else:
            iso = ",".join([str(x) for x in isotopologues])
        options = {"iso_ids_list" : iso, "numin" : lower_bound, "numax" : upper_bound}
        opts = "&".join(["{}={}".format(key, value) for key, value in options.items()])
        url = "{}{}".format(url, opts)
        parameters = []
        id = self.molecule_ids[molecule]
        for record in self.molecular_line_records(urlopen(url)):
            line = HitranSpectralLine()
            line.parse_ascii_record(record, id, lower_bound, upper_bound)
            parameters.append(line)
        return parameters

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
        parameters = []
        for record in cursor.fetchall():
            line = HitranSpectralLine()
            for attr, data in zip(members.values(), record):
                setattr(line, attr, data)
            parameters.append(line)
        connection.close()

    @staticmethod
    def molecular_line_records(request):
        """Reads the next record from the html molecular line parameter table.

        Args:
            request: A urllib.request.Request object.

        Yields:
            Record of the HITRAN database.
        """
        record_byte_length = 162
        while True:
            record = request.read(record_byte_length).decode("utf-8")
            if record == "":
                return
            yield record

    def read_from_ascii(self, path, molecule, isotopologue_ids="all", lower_bound=0.,
                        upper_bound=finfo("f4").max):
        """Reads HITRAN molecular line parameters from an ascii file.

        Args:
            path: Path to ascii file.
            molecule: String molecule chemical formula (e.g., "H2O").
            isotopologues: List of HITRAN isotopologue ids.
            lower_bound: Lower bound of spectral range [cm-1], inclusive.
            upper_bound: Upper bound of spectral range [cm-1], inclusive.
        """
        molecule_id = self.molecule_ids[molecule]
        if isotopologue_ids == "all":
            isotopologue_ids = [x+1 for x in range(len(self.isotopologues[molecule]))]
        else:
            isotopologue_ids = [int(x) for x in isotopologues]
        parameters = []
        id = self.molecule_ids[molecule]
        with open(path, "r") as file:
            for record in file:
                try:
                    line = HitranSpectralLine()
                    line.parse_ascii_record(record, id, lower_bound, upper_bound)
                    parameters.append(line)
                except RecordError:
                    continue
        return parameters


def create_database(database, molecule, parameters):
    """Creates/ingests data into a SQLite database.

    Args:
        database: Path to SQLite database that will be create/added to.
        molecule: Chemical formula of molecule.
        parameters: List of HitranSpectralLine objects.
    """
    connection = connect(database)
    cursor = connection.cursor()
    name = scrub(molecule)
    types = {int : "INTEGER", float32 : "REAL"}
    headings = ", ".join(["{} {}".format(x.name, types[x.type]) for x in columns[1:]])
    cursor.execute("CREATE TABLE {} ({})".format(name, headings))
    value_subst = ", ".join(["?" for _ in columns[1:]])
    types = {int: int, float32 : float64}
    for p in parameters:
        values = tuple([types[x.type](getattr(p, members[x.name])) for x in columns[1:]])
        cursor.execute("INSERT INTO {} VALUES ({})".format(name, value_subst), values)
    connection.commit()
    connection.close()


def write_to_ascii(path, id, parameters):
    """Writes HITRAN line parameters to an ascii file.

    Args:
        path: Path to output file.
        id: HITRAN identifier for the molecule.
        parameters: List of HitranSpectralLine objects.
    """
    formats = {int : "d", float32 : "g"}
    with open(path, "w") as file:
        for line in parameters:
            record = []
            for column in columns:
                s = "{{:>{}{}}}".format(column.width, formats[column.type])
                if column.name == "hitran_id":
                    record.append(s.format(id))
                else:
                    datum = getattr(line, members[column.name])
                    if column.name == "isotopologue" and datum > 9:
                        if datum == 10:
                            value = "0"
                        else:
                            value = chr(ord("A") + datum - 11)
                    else:
                        value = s.format(datum)
                    record.append(value[len(value)-column.width:])
                    if column.name == "strength":
                        record.append(" "*10)
            file.write("{}\n".format("".join(record)))
