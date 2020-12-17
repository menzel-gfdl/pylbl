from logging import getLogger

from netCDF4 import Dataset
from numpy import copy, exp, identity, int32, log, nonzero, ones, power, sqrt, where, zeros
from scipy.special import factorial

from .line_mixing_fortran import calculate_coefficients, create_relaxation_matrix
from .line_parameters import c2, reference_temperature
from ..tips import TIPS_REFERENCE_TEMPERATURE
from .wigner3j import wigner_3j


info = getLogger(__name__).info


class LineMixing(object):
    def __init__(self, path):
        """Reads the (w0, b0) fitted parameters taken from doi: 10.1016/j.jqsrt.2014.09.017.

        Args:
            path - Path to netcdf file.
        """
        info("Reading line mixing data from {}.".format(path))
        with Dataset(path, "r") as dataset:
            for branch in ["pp", "pq", "pr", "qp", "qq", "qr", "rp", "rq", "rr"]:
                setattr(self, "w0{}".format(branch),
                        copy(dataset.variables["w0{}".format(branch)], order="F"))
                setattr(self, "b0{}".format(branch),
                        copy(dataset.variables["b0{}".format(branch)], order="F"))

    def create_relaxation_matrix(self, li, lf, ji, jf, population, temperature, iso,
                                 gamma, dk0, mask, I, I_off):
        """Calculates the relaxation matrix as in doi: 10.1016/j.jqsrt.2014.09.017.

        Args:
            li:
            lf:
            ji:
            jf:
            population:
            temperature: Temperature [K].
            iso: Isotopologue id.
            gamma:
            dk0:
            mask:
            I: Identity matrix.
            I_off: Matrix with zeros on-diagonal and ones off-diagonal.

        Returns:
            Relaxation matrix.
        """
        tlog = log(reference_temperature/temperature)
        w = zeros((ji.size, ji.size))
        if li <= lf:
            a, b, j1, j2 = li, lf, ji[:], jf[:]
        else:
            a, b, j1, j2 = lf, li, jf[:], ji[:]
        #xmask = where((j1.reshape(j1.size, 1) > j1), 0., 1.)*mask
        xmask = where((j1 > j1.reshape(j1.size, 1)), 0., 1.)*mask
        for i in range(ji.size):
            branch1 = branch(j1[i], j2[i])
            w0, b0 = zeros(ji.size), zeros(ji.size)
            indices = nonzero(xmask[i, :])[0]
            for j in indices:
                branch2 = branch(j1[j], j2[j])
                w0[j] = getattr(self, "w0{}{}".format(branch1, branch2))[a, b, j1[i], j1[j]]
                b0[j] = getattr(self, "b0{}{}".format(branch1, branch2))[a, b, j1[i], j1[j]]
            ycal = exp(w0[indices] - b0[indices]*tlog)
            w[indices, i] = ycal[:]
            w[i, indices] = ycal[:]*population[i]/population[indices]

        # Set diagonal terms to the line-broadening parameters.
        w = -1.*abs(w)*I_off + gamma*I

        for i in range(gamma.size):
            sum_lw = sum(mask[i, i+1:]*abs(dk0[i+1:])*w[i+1:, i])
            sum_up = sum(mask[i, :i+1]*abs(dk0[:i+1])*w[:i+1, i])
            if sum_lw == 0.:
                w[i+1:, i] = 0.
                w[i, i+1:] = 0.
            else:
                w[i+1:, i] = -1.*w[i+1:, i]*sum_up/sum_lw
                w[i, i+1:] = w[i+1:, i]*population[i]/population[i+1:]
        return w

    def first_order_coefficients(self, spectral_lines, temperature, pressure, mixing_ratio):
        """Calculates 1st order line-mixing coefficients.

        Args:
            spectral_lines: SpectralLines object.
            temperature: Temperature [K].
            pressure: Pressure [atm].
            mixing_ratio: Mixing ratio.
        """
        bands = create_bands(spectral_lines)
        y = zeros(spectral_lines.v.size)
        for i, indices in enumerate(bands.values()):
            info("Calculating line-mixing coefficents (band {}, {} lines)".format(i, len(indices)))

            #Determine the initial and final l values.
            li = int32(spectral_lines.q2[indices[0]]["l2"])
            lf = int32(spectral_lines.q1[indices[0]]["l2"])
            if li > 8 or lf > 8 or abs(li - lf) > 1:
                continue

            #Create array of initial and final j values.
            ji, jf = zeros(len(indices), dtype=int), zeros(len(indices), dtype=int32)
            for j, (qns2, qns1) in enumerate(zip(spectral_lines.q2[indices],
                                                 spectral_lines.q1[indices])):
                ji[j] = qns2["j"]
                jf[j] = qns1["j"]

            #Adjust the level populations for temperature.
            n = temperature_adjust_population(spectral_lines.iso[indices[0]],
                                              spectral_lines.population[indices],
                                              spectral_lines.en[indices],
                                              temperature,
                                              spectral_lines.q.total_partition_function)

            #Adjust the doppler halfwidth for temperature.
            gamma = temperature_adjust_halfwidth(spectral_lines.n_air[indices],
                                                 spectral_lines.n_self[indices],
                                                 spectral_lines.gamma_air[indices],
                                                 spectral_lines.gamma_self[indices],
                                                 mixing_ratio,
                                                 temperature)

            #Sort parameters in decreasing order of intensity.
            s = spectral_lines.v[indices]*n[:]*power(spectral_lines.dipole[indices], 2)
            sorted_indices = s.argsort()
            ji_s = ji[sorted_indices[::-1]]
            jf_s = jf[sorted_indices[::-1]]
            s_s = s[sorted_indices[::-1]]
            v_s = spectral_lines.v[indices][sorted_indices[::-1]]
            gamma_s = gamma[sorted_indices[::-1]]
            dair_s = spectral_lines.d_air[indices][sorted_indices[::-1]]
            n_s = n[sorted_indices[::-1]]
            dk0_s = spectral_lines.dk0[indices][sorted_indices[::-1]]
            dipole_s = spectral_lines.dipole[indices][sorted_indices[::-1]]

            #Create some useful constants.
            I = identity(gamma_s.size)
            unity = ones((gamma_s.size, gamma_s.size))
            I_off = unity - I
            x = mask(spectral_lines.iso[indices[0]], ji_s)

            #Create relaxation matrix.
#           w = self.create_relaxation_matrix(li, lf, ji_s, jf_s, n_s, temperature,
#                                             spectral_lines.iso[indices[0]], gamma_s,
#                                             dk0_s, x, I, I_off)
            w = zeros((ji_s.size, ji_s.size), order="F")
            create_relaxation_matrix(ji_s.size, temperature, spectral_lines.iso[indices[0]],
                                     li, lf, gamma_s, n_s, dk0_s, ji_s, jf_s,
                                     self.b0pp, self.b0pq, self.b0pr,
                                     self.b0qp, self.b0qq, self.b0qr,
                                     self.b0rp, self.b0rq, self.b0rr,
                                     self.w0pp, self.w0pq, self.w0pr,
                                     self.w0qp, self.w0qq, self.w0qr,
                                     self.w0rp, self.w0rq, self.w0rr, w)

            #Calculate first order line-mixing coefficients.
#           y[indices] = calculate_coefficients(w, dipole_s, v_s, I_off, x) * \
#                        pressure
            yt = zeros(ji_s.size)
            calculate_coefficients(ji_s.size, spectral_lines.iso[indices[0]],
                                   dipole_s, ji_s, v_s, w, yt)
            y[indices] = yt[:]*pressure

        return y


def temperature_adjust_population(iso, population, en, temperature, partition_function):
    """Adjusts level populations to input temperature.

    Args:
        iso: Isotopologue id.
        population: Upper level population.
        en: Lower state enegery [cm-1].
        temperature: Temperature [K].
        partition_function: Total partition function routine.

    Returns:
        Temperature-corrected population.
    """
    return population*partition_function(TIPS_REFERENCE_TEMPERATURE, iso) * \
           exp(c2*en*(1./temperature - 1./TIPS_REFERENCE_TEMPERATURE)) / \
           partition_function(temperature, iso)


def temperature_adjust_halfwidth(n_air, n_self, gamma_air, gamma_self, mixing_ratio,
                                 temperature):
    """Adjusts input halfwidths to input temperature.

    Args:
        n_air: Air-broadened halfwidth temperature dependence power.
        n_self: Self-broadened halfwidth temperature dependence power.
        gamma_air: Air-broadened halfwidth at reference temperature [cm-1].
        gamma_self: Self-broadened halfwidth at reference temperature [cm-1].
        mixing_ratio: Molecular mixing ratio.
        temperature: Temperature [K].

    Returns:
        Temperature-corrected halfwidth [cm-1].
    """
    return (1. - mixing_ratio)*gamma_air*power(reference_temperature/temperature, n_air) + \
           mixing_ratio*gamma_self*power(reference_temperature/temperature, n_self)


def create_bands(spectral_lines):
    """Groups the lines in the same vibrational bands.

    Args:
        spectral_lines: SpectralLines object.

    Returns:
        Dictionary of lists of indices.
    """
    bands = {}
    for i, (upper, lower) in enumerate(zip(spectral_lines.q1, spectral_lines.q2)):
        key = tuple([upper[x] for x in ["v1", "v2", "l2", "v3"]] +
                    [lower[x] for x in ["v1", "v2", "l2", "v3"]])
        if key in bands:
            bands[key].append(i)
        else:
            bands[key] = [i]
    return bands


def branch(ji, jf):
    """Returns the branch designation string for the transition.

    Args:
        ji:
        jf:

    Returns:
        Character describing transition.
    """
    if ji > jf:
        return "p"
    elif ji == jf:
        return "q"
    else:
        return "r"


def mask(iso, ji):
    """Creates a mask for transitions.

    Args:
        iso: Isotopologue id.
        ji:

    Returns:
        Array of zeros and ones used to mask out transitions.
    """
    if iso > 2 and iso != 7 and iso != 10:
        return where((ji.reshape((ji.size, 1)) - ji) % 2 != 0, 0., 1.)
    else:
        return ones((ji.size, ji.size))


def rigid_rotor_dipole_matrix_element(ji, jf, l2i, l2f):
    """Computes a rigid rotor dipole matrix element from equation 3 in
       doi: 10.1016/j.jqsrt.2004.11.011

    Args:
        ji: Initial state J value.
        jf: Final state J value.
        l2i: Initial state l2 value.
        l2f: Final state l2 value.

    Returns:
        Rigid rotor dipole matrix element.
    """
    return power(-1., jf + l2f + 1)*sqrt(2*jf + 1) * \
        wigner_3j(ji, 1, jf, l2i, l2f - l2i, -l2f)


#def calculate_coefficients(relaxation_matrix, dipole, line_center, I_off, mask):
#    """Calculates first-order line-mixing coefficients using equation 6 of
#       doi: 10.1016/j.jqsrt.2004.11.011.
#
#    Args:
#        relaxation_matrix: Relaxation matrix.
#        line_center: Line center wavenumber [cm-1].
#        I_off: 2D matrix with 0s on diagnal and 1s off-diagonal.
#        mask: Matrix used to mask out certain transitions.
#
#    Returns:
#        First-order line-mixing coefficients.
#    """
#    y = zeros(line_center.size)
#    dipole_ratio = dipole/dipole.reshape((dipole.size, 1))
#    x = line_center.reshape((line_center.size, 1)) - line_center
#    dw = where(abs(x) < 1.e-4, 1.e-4, x)
#    x = mask*I_off*dipole_ratio/dw
#    for i in range(line_center.size):
#        y[i] = 2.*sum(x[i, :]*relaxation_matrix[:, i])
#    return y
