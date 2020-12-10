from collections import namedtuple, OrderedDict
from re import match

from numpy import float64


def linear_molecule_quantum_numbers(value):
    m = match(r"ElecStateLabel=X;v1=([0-9]+);v2=([0-9]+);l2=([0-9]+);v3=([0-9]+);J=([0-9]+);",
              value)
    if not m:
        raise ValueError("invalid quantum numbers in {}".format(value))
    return OrderedDict((x, int(y)) for x, y in zip(["v1", "v2", "l2", "v3", "j"],
                                                   m.groups()))


HitranParameter = namedtuple("HitranParameter", ["api_name", "longname", "shortname",
                                                 "dtype"])

_species = ["air", "co2", "he", "h2", "self"]

_base_parameters = {
    "center": HitranParameter(api_name="nu",
                              longname="center",
                              shortname="v",
                              dtype=float64),
    "elower": HitranParameter(api_name="elower",
                              longname="lower_state_energy",
                              shortname="en",
                              dtype=float64),
    "id": HitranParameter(api_name="molec_id",
                          longname="hitran_id",
                          shortname="id",
                          dtype=int),
    "iso": HitranParameter(api_name="local_iso_id",
                           longname="isotopologue",
                           shortname="iso",
                           dtype=int),
    "strength": HitranParameter(api_name="sw",
                                longname="strength",
                                shortname="s",
                                dtype=float64),
    "trans_id": HitranParameter(api_name="trans_id",
                                longname="transition_id",
                                shortname="line_id",
                                dtype=int),
    "upper_level_quantum_numbers": HitranParameter(api_name="statep",
                                                   longname="upper_level_quantum_numbers",
                                                   shortname="q1",
                                                   dtype=linear_molecule_quantum_numbers),
    "lower_level_quantum_numbers": HitranParameter(api_name="statepp",
                                                   longname="lower_level_quantum_numbers",
                                                   shortname="q2",
                                                   dtype=linear_molecule_quantum_numbers),
    "lower_level_degeneracy": HitranParameter(api_name="gpp",
                                              longname="lower_level_degeneracy",
                                              shortname="g",
                                              dtype=int)
}

_beta_parameters = {
    "beta_g_{}".format(x): HitranParameter(api_name="beta_g_{}".format(x),
                                           longname="galatry profile coefficients for {}-broadening".format(x),
                                           shortname="beta_g_{}".format(x),
                                           dtype=float64) for x in ["air", "self"]
}

_delta_parameters = {
    "delta_{}".format(x): HitranParameter(api_name="delta_{}".format(x),
                                          longname="{}_broadened_pressure_shift".format(x),
                                          shortname="d_{}".format(x),
                                          dtype=float64) for x in _species
}

_deltap_parameters = {
    "deltap_{}".format(x): HitranParameter(api_name="deltap_{}",
                                           longname="",
                                           shortname="deltap_{}".format(x),
                                           dtype=float64) for x in ["air", "h2", "self"]
}

_gamma_parameters = {
    "gamma_{}".format(x): HitranParameter(api_name="gamma_{}".format(x),
                                          longname="{}_broadened_halfwidth".format(x),
                                          shortname="gamma_{}".format(x),
                                          dtype=float64) for x in _species
}

_n_parameters = {
    "n_{}".format(x): HitranParameter(api_name="n_{}".format(x),
                                      longname="{}_broadened_temperature_dependence".format(x),
                                      shortname="n_{}".format(x),
                                      dtype=float64) for x in _species
}

PARAMETERS = {**_base_parameters, **_beta_parameters, **_delta_parameters,
              **_deltap_parameters, **_gamma_parameters, **_n_parameters}
