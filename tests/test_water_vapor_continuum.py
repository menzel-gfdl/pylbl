from numpy import arange

from pyrad.lbl.continua import WaterVaporContinuum


if __name__ == "__main__":
    water = WaterVaporContinuum()
    grid = arange(1., 3250., 0.01)
    print(water.absorption_coefficient(300., 101300., 0.02, grid))
