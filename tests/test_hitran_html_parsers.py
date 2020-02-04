from unittest import main, TestCase
from urllib.request import urlopen

from pyrad.lbl.hitran_html_parsers import HitranIsotopologueHTMLParser, HitranMoleculeIdHTMLParser


class TestHitranHTMLParsers(TestCase):

    def test_molecule_id(self):
        ids = ("H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2",
               "NH3", "HNO3", "OH", "HF", "HCl", "HBr", "HI", "ClO", "OCS", "H2CO",
               "HOCl", "N2", "HCN", "CH3Cl", "H2O2", "C2H2", "C2H6", "PH3", "COF2",
               "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2", "NO+", "HOBr", "C2H4",
               "CH3OH", "CH3Br", "CH3CN", "CF4", "C4H2", "HC3N", "H2", "CS", "SO3",
               "C2N2", "COCl2", "SO", "C3H4", "CH3", "CS2")
        parser = HitranMoleculeIdHTMLParser()
        parser.feed(urlopen("https://hitran.org/docs/molec-meta/").read().decode("utf-8"))
        for i, molecule in enumerate(ids):
            if molecule in ["SO", "C3H4", "CH3", "CS2"]:
                continue
            self.assertEqual(i+1, parser.molecule_ids[molecule])

    def test_isotopologue_id(self):
        molecule = "O3"
        isotopologue_ids = [16, 17, 18, 19, 20]
        parser = HitranIsotopologueHTMLParser()
        parser.feed(urlopen("https://hitran.org/docs/iso-meta/").read().decode("utf-8"))
        for i, iso in enumerate(parser.isotopologues[molecule]):
            self.assertEqual(iso.id, isotopologue_ids[i])


if __name__ == "__main__":
    main()
