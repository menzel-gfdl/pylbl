from unittest import main, TestCase

from pyrad.lbl.hitran.molecules import molecules


class TestMolecules(TestCase):

    def test_molecules(self):
        formulae = ("H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2",
                    "NH3", "HNO3", "OH", "HF", "HCl", "HBr", "HI", "ClO", "OCS", "H2CO",
                    "HOCl", "N2", "HCN", "CH3Cl", "H2O2", "C2H2", "C2H6", "PH3", "COF2",
                    "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2", "NO+", "HOBr", "C2H4",
                    "CH3OH", "CH3Br", "CH3CN", "CF4", "C4H2", "HC3N", "H2", "CS", "SO3",
                    "C2N2", "COCl2", "SO", "C3H4", "CH3", "CS2")
        ids = molecules("https://hitran.org/docs/molec-meta/")
        for i, formula in enumerate(formulae):
            if formula not in ["SO", "C3H4", "CH3", "CS2"]:
                self.assertEqual(i+1, ids[formula])


if __name__ == "__main__":
    main()
