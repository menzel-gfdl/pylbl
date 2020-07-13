from os.path import abspath
from tempfile import NamedTemporaryFile
from unittest import main, TestCase
from urllib.parse import quote
from urllib.request import urlopen

from pyrad.utils.database_utilities import ascii_table_records, scrub


class TestDatabaseUtilities(TestCase):

    def test_ascii_table_records(self):
        test_data = ["this is line {}.".format(x) for x in range(1000)]
        with NamedTemporaryFile() as test_file:
            test_file.write("\n".join(test_data).encode("utf-8"))
            url = urlopen("file://" + quote(abspath(test_file.name)))
            for i, record in enumerate(ascii_table_records(url)):
                self.assertEqual(record, test_data[i])

    def test_scrub(self):
        self.assertEqual(scrub("foo; DROP TABLE bar;"), "foo")


if __name__ == "__main__":
    main()
