from copy import copy as shallow_copy

from numpy import asarray, copy, exp, searchsorted, sqrt, zeros

from .line_mixing import rigid_rotor_dipole_matrix_element
from .line_parameters import c2, reference_temperature
from ..tips import TIPS_REFERENCE_TEMPERATURE


class SpectralLines(object):
    """Hitran spectral line parameters for a molecule.

    Attributes:
        delta: Numpy array of air-broadended pressure shifts [cm-1 atm-1] (lines).
        dipole: Numpy array of transition moments [cm] (lines).
        dk0: Numpy array of rigid rotor dipole matrix elements (lines).
        en: Numpy array of transition lower state energies [cm-1] (lines).
        g: Numpy array of lower state statistical weights (lines).
        gamma_air: Numpy array of air-broadened halfwidths [cm-1 atm-1] (lines).
        gamma_self: Numpy array of self-broadened halfwidths [cm-1 atm-1] (lines).
        id: HITRAN molecule id.
        iso: Numpy array of HITRAN isotopologue ids (lines).
        mass: Numpy array of isotopologue masses [g] (lines).
        n: Numpy array of air-broadened temperature dependence powers (lines).
        population: Numpy array of lower state populations (lines).
        s: Numpy array of line strengths [cm] (lines).
        q: TotalPartitionFunction object.
        v: Numpy array of transition wavenumbers [cm-1] (lines).
    """

    def __init__(self, database, total_partition_function):
        """Sorts line parameters by transition wavenumber, and partially corrects the line strengths.

        Args:
            database: Hitran object.
            total_partition_function: TotalPartitionFunction object.

        Raises:
            EmptySpectraError: No molecular line parameters are detected.
        """
        # Create member arrays.
        for x in database.parameters:
            setattr(self, x.shortname, copy(getattr(database, x.shortname)))

        # Correct for Hitran counting weirdness (1, 2, 3, ... 9, 0, a, b, ...)
        for i in range(self.iso.size):
            if self.iso[i] == 0:
                self.iso[i] = 10

        # Get the mass of the isotopologues.
        self.mass = asarray([float(database.isotopologues[x-1].mass) for x in self.iso])
        self.line_profile = database.line_profile
        self.q = total_partition_function

        # Partially correct line strengths.
        self.s[:] *= self.temperature_correct_line_strength(self.q, TIPS_REFERENCE_TEMPERATURE,
                                                            self.iso, self.en, self.v)

        if self.line_profile.line_mixing is not None:
            self.population = self.g[:]*exp(c2*self.en[:]/reference_temperature) / \
                              self.q.total_partition_function(TIPS_REFERENCE_TEMPERATURE, self.iso[:])
            self.dipole = sqrt(self.s[:]/(self.g[:]*self.v[:]))
            self.dk0 = zeros(self.q2.size)
            for i in range(self.q2.size):
                self.dk0[i] = rigid_rotor_dipole_matrix_element(self.q2[i]["j"],
                                                                self.q1[i]["j"],
                                                                self.q2[i]["l2"],
                                                                self.q1[i]["l2"])

    def absorption_coefficient(self, temperature, pressure, partial_pressure, wavenumber,
                               cut_off=25.):
        """Calculates the absorption coefficient.

        Args:
            temperature: Temperature [K].
            pressure: Pressure [atm].
            partial_pressure: Partial pressure [atm].
            wavenumber: Numpy array of wavenumbers [cm-1] (wavenumber).
            cut_off: Distance [cm-1] from the transition frequency where the line is cut off.

        Returns:
            Numpy array of absorption coefficients [cm2] (wavenumber).
        """
        lines = shallow_copy(self)
        profile = shallow_copy(lines.line_profile)
        profile.update(lines, temperature, pressure, partial_pressure)
        lines.s = lines.correct_line_strengths(temperature)
        lines.v = lines.pressure_shift_transition_wavenumbers(pressure)
        k = zeros(wavenumber.size)
        for i in range(lines.s.size):
            left = searchsorted(wavenumber, lines.v[i] - cut_off, side="left")
            right = searchsorted(wavenumber, lines.v[i] + cut_off, side="right")
            k[left:right] += lines.s[i]*profile.profile(lines, wavenumber[left:right], i)
        return k

    def correct_line_strengths(self, temperature):
        """Temperature-corrects the line strengths.

        Args:
            temperature: Temperature [K].

        Returns:
            Numpy array of corrected line strengths [cm] (lines).
        """
        return self.s[:]*(1./self.temperature_correct_line_strength(self.q, temperature,
                                                                    self.iso, self.en,
                                                                    self.v))

    def pressure_shift_transition_wavenumbers(self, pressure):
        """Pressure-shifts transition wavenumbers.

        Args:
            pressure: Pressure [atm].

        Returns:
            Numpy array of pressure-shifted transition wavenumbers [cm-1] (lines).
        """
        return self.v[:] + self.d_air[:]*pressure

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
            Temperature correction factor.
        """
        # Divide-by-zeros may occur for transition wavenumbers close to zero, like those
        # for the O16-O17 isotopologue of O2.
        return q.total_partition_function(t, iso[:])/(exp(c2*en[:]/t)*(1. - exp(c2*v[:]/t)))
