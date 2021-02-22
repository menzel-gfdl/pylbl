import matplotlib.pyplot as plt
from numpy import arange, seterr

from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA


if __name__ == "__main__":
    seterr(all="raise")
    grid = arange(1., 3250., 0.01)
    for interaction in [["CO2", "CO2"], ["N2", "N2"], ["O2", "N2"], ["O2", "O2"]]:
        mol = HitranCIA(*interaction)
        for temperature in [211.11, 235.45, 250.23, 270.45, 290.11, 300.1]:
            xsec = mol.absorption_coefficient(temperature, grid)
            plt.plot(grid, xsec, label="{}".format(temperature))
        plt.title("-".join(interaction))
        plt.legend()
