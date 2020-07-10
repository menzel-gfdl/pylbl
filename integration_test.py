from logging import basicConfig, INFO

import matplotlib.pyplot as plt
from numpy import seterr

from pyrad.gas_optics import Gas
from pyrad.lbl.hitran.database import Database
from pyrad.lbl.hitran.lines import SpectralLines
from pyrad.lbl.line_profiles import Voigt
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.utils.grids import UniformGrid1D


if __name__ == "__main__":
    basicConfig(format="%(asctime)-15s - %(funcName)s(%(lineno)d): %(message)s", level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")

    #Configure layer.
    mb_to_atm = 0.000986923
    spectral_grid = UniformGrid1D(1., 500., 0.1)
    temperature = 288.99 #K
    pressure = 983.88*mb_to_atm #atm
    water_vapor_abundance = 0.006637074
    partial_pressure = pressure*water_vapor_abundance #atm

    #Calculate absorption coefficient.
    hitran = Database()
    water = Gas("H2O", hitran)
    k = water.absorption_coefficient(temperature, pressure, partial_pressure,
                                     spectral_grid.points, Voigt())

    #Plot absorption coefficient.
    plt.plot(spectral_grid.points, k)
    plt.show()
