from logging import info

from numpy import argsort, asarray, exp, log, power, searchsorted, sqrt, sum, zeros

from .database import PARAMETERS
from ..line_profiles import Doppler, Lorentz, Voigt
from ..tips import TIPS_REFERENCE_TEMPERATURE


class _SpectralLine(object):
    """HITRAN spectral parameters for a single molecular line.

    Attributes:
        delta: Air-broadended pressure shift [cm-1 atm-1].
        en: Transition lower state energy [cm-1].
        gamma_air: Air-broadened halfwidth [cm-1 atm-1].
        gamma_self: Self-broadened halfwidth [cm-1 atm-1].
        id: HITRAN molecule id.
        iso: HITRAN isotopologue id.
        line_id: HITRAN transition id.
        n: Air-broadened temperature dependence power.
        s: Line strength [cm-1].
        v: Transition wavenumber [cm-1].
    """
    def __init__(self):
        for parameter in PARAMETERS:
            setattr(self, parameter.shortname, None)

    def parse_ascii_record(self, record):
        """Parses spectral line parameters from a HITRAN database record.

        Args:
            record: HITRAN database record.

        Returns:
            The current _SpectralLine object.
        """
        data = [x.strip() for x in record.split(",")]
        for i, parameter in enumerate(PARAMETERS):
            setattr(self, parameter.shortname, parameter.type(data[i]))
        if self.iso == 0:
            self.iso = 10
        return self


class SpectralLines(object):
    """Hitran spectral line parameters for a molecule.

    Attributes:
        delta: Numpy array of air-broadended pressure shifts [cm-1 atm-1] (lines).
        en: Numpy array of transition lower state energies [cm-1] (lines).
        gamma_air: Numpy array of air-broadened halfwidths [cm-1 atm-1] (lines).
        gamma_self: Numpy array of self-broadened halfwidths [cm-1 atm-1] (lines).
        id: HITRAN molecule id.
        iso: Numpy array of HITRAN isotopologue ids (lines).
        mass: Numpy array of isotopologue masses (lines).
        n: Numpy array of air-broadened temperature dependence powers (lines).
        s: Numpy array of line strengths [cm-1] (lines).
        q: TotalPartitionFunction object.
        v: Numpy array of transition wavenumbers [cm-1] (lines).
    """

    def __init__(self, molecule, isotopologues, database, total_partition_function):
        """Sorts line parameters by transition wavenumber, and partially corrects the line strengths.

        Args:
            molecule: String chemical formula.
            isotopologues: List of isotopologue namedtuples.
            database: Database object.
            total_partition_function: TotalPartitionFunction object.

        Raises:
            EmptySpectraError: No molecular line parameters are detected.
        """
        #Download spectral line parameters and sort them by transition wavenumber.
        lines = []
        for record in database.records(molecule, isotopologues):
            try:
                lines.append(_SpectralLine().parse_ascii_record(record))
            except ValueError as error:
                if not "could not convert string" in str(error):
                    raise
        if not lines:
            raise EmptySpectraError("no spectral lines in the input spectral range.")
        parameters = sorted(lines, key=lambda line: line.v)
        num_lines = len(lines)
        info("Found {} spectral lines.".format(num_lines))

        #Create member arrays.
        for key, value in vars(parameters[0]).items():
            setattr(self, key, zeros(num_lines, dtype=type(value)))

        #Convert from array-of-structs to struct-of-arrays.
        for i, p in enumerate(parameters):
            for key, value in vars(p).items():
                self.__dict__[key][i] = value

        #Partially correct line strengths.
        self.q = total_partition_function
        self.s[:] *= self.temperature_correct_line_strength(self.q, TIPS_REFERENCE_TEMPERATURE,
                                                            self.iso[:], self.en[:], self.v[:])

        #Get the mass of the isotopologues.
        self.mass = asarray([database.isotopologues[molecule][x-1].mass for x in self.iso])

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
