from unittest import main, TestCase

from pylbl.database_utilities import scrub


class TestTips(TestCase):

    def test_scrub(self):
        self.assertEqual(scrub("foo; DROP TABLE bar;"), "foo")


if __name__ == "__main__":
    main()
