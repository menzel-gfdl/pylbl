from unittest import main, TestCase

from pyrad.lbl.hitran.database import Database


class TestDatabase(TestCase):

    def test_database(self):
        database = Database()
        formulae = ["H2O", "CO2", "O3", "N2O", "CH4"]
        for formula in formulae:
            for record in database.records(formula):
                pass


if __name__ == "__main__":
    main()
