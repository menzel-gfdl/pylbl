from collections import namedtuple
from logging import info
from urllib.request import urlopen

from numpy import float32

from .isotopologues import isotopologues
from .molecules import molecules
from ...utils.database_utilities import ascii_table_records


_Parameter = namedtuple("_Parameter", ["api_name", "longname", "shortname", "type"])


PARAMETERS = [
    _Parameter(api_name="molec_id", longname="hitran_id", shortname="id", type=int),
    _Parameter(api_name="local_iso_id", longname="isotopologue", shortname="iso", type=int),
    _Parameter(api_name="nu", longname="center", shortname="v", type=float32),
    _Parameter(api_name="sw", longname="strength", shortname="s", type=float32),
    _Parameter(api_name="gamma_air", longname="air_broadened_halfwidth",
               shortname="gamma_air", type=float32),
    _Parameter(api_name="gamma_self", longname="self_broadened_halfwidth",
               shortname="gamma_self", type=float32),
    _Parameter(api_name="elower", longname="lower_state_energy", shortname="en",
               type=float32),
    _Parameter(api_name="n_air", longname="air_broadened_temperature_dependence",
               shortname="n", type=float32),
    _Parameter(api_name="delta_air", longname="air_broadened_pressure_shift",
               shortname="delta", type=float32),
    _Parameter(api_name="trans_id", longname="transition_id", shortname="line_id", type=int)]


class Database(object):
    """HITRAN database front-end.

    Attributes:
       molecule_ids: Dictionary mapping chemical formulae to HITRAN molecule ids.
       isotopologues: Dictionary mapping chemical formuale to lists of Isotopologue objects.
    """

    def __init__(self):
        """Creates dictionaries of HITRAN ids by scraping the HITRAN website."""
        self.molecule_ids = molecules("https://hitran.org/docs/molec-meta/")
        self.isotopologues = isotopologues("https://hitran.org/docs/iso-meta/")

    def records(self, molecule, isotopologues="all", lower_bound=0.,
                upper_bound=10.e6, parameters=PARAMETERS):
        """Downloads HITRAN molecular line parameters from the web.

        Args:
            molecule: String molecule chemical formula (e.g., "H2O").
            isotopologues: List of HITRAN isotopologue ids.
            lower_bound: Lower bound of spectral range [cm-1], inclusive.
            upper_bound: Upper bound of spectral range [cm-1], inclusive.
        """
        if isotopologues == "all":
            isotopologues = self.isotopologues[molecule]
        options = {"iso_ids_list" : ",".join([str(x.id) for x in isotopologues]),
                   "numin" : lower_bound, 
                   "numax" : upper_bound,
                   "request_params" : ",".join([x.api_name for x in parameters])}
        html_options = "&".join(["{}={}".format(key, value) for key, value in options.items()])
        url = "http://hitran.org/lbl/api?{}".format(html_options)
        info("Downloading Hitran database for {} from {}.".format(molecule, url))
        for record in ascii_table_records(urlopen(url)):
            yield record
