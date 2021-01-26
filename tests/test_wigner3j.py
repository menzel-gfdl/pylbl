from unittest import main, TestCase

from numpy import seterr, sqrt

from pyrad.lbl.hitran.wigner3j import wigner_3j


class TestWigner3j(TestCase):

    def test_wigner_3j(self):
        inputs = [((24, 1, 23, 0, 1, -1), (1/7)*sqrt(23/94)),  # Special case 3.
                  ((4, 1, 3, 1, 0, -1), -0.5*sqrt(5/21)),
                  ((18, 1, 19, 2, 0, -2), -sqrt(119/9139)),
                  ((23, 1, 22, 3, -1, -2), (1/3)*sqrt(65/1081)),
                  ((39, 1, 38, 3, 1, -4), sqrt(30/11297)),
                  ((1, 1, 1, 1, 0, -1), -sqrt(1/6)),  # Special case 4.
                  ((1, 1, 1, 1, -1, 0), sqrt(1/6)),
                  ((2, 1, 2, 2, 0, -2), -sqrt(2/15)),
                  ((3, 1, 3, 3, 0, -3), -0.5*sqrt(3/7)),
                  ((4, 1, 4, 4, 0, -4), (-2/3)*sqrt(0.2)),
                  ((5, 1, 5, 5, 0, -5), -sqrt(5/66)),
                  ((2, 1, 2, 0, 1, -1), -sqrt(0.1)),  # General case.
                  ((10, 1, 10, 0, 1, -1), -sqrt(1/42)),
                  ((24, 1, 24, 0, 1, -1), -(1/7)*sqrt(0.5)),
                  ((2, 1, 2, 1, 0, -1), sqrt(1/30)),
                  ((22, 1, 22, 1, 0, -1), (1/3)*sqrt(1/2530)),
                  ((12, 1, 12, 1, -1, 0), -0.2*sqrt(0.5)),
                  ((7, 1, 7, 3, -1, -2), 0.5*sqrt(5/42)),
                  ((29, 1, 29, 3, -1, -2), 6*sqrt(2/8555))]
        for x in inputs:
            self.assertAlmostEqual(wigner_3j(*x[0]), x[1], places=7)


if __name__ == "__main__":
    seterr(divide="raise", over="raise", invalid="raise")
    main()
