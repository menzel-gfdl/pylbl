from unittest import main, TestCase

from pyrad.lbl.hitran.isotopologues import isotopologues


class TestIsotopologues(TestCase):

    def test_isotopologues(self):
        formula = "O3"
        ids = [16, 17, 18, 19, 20]
        iso = isotopologues("https://hitran.org/docs/iso-meta/")
        for i, id_ in enumerate(ids):
            self.assertEqual(id_, iso[formula][i].id)


if __name__ == "__main__":
    main()
