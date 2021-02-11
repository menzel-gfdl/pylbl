from numpy import arange

from pyrad.lbl.continua import OzoneContinuum


if __name__ == "__main__":
    continuum = OzoneContinuum()
    grid = arange(9000., 50000., 0.1)
    print(continuum.absorption_coefficient(grid))
