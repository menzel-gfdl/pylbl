import matplotlib.pyplot as plt
from numpy import arange, seterr

from pyrad.lbl.hitran.cross_sections import HitranCrossSection


if __name__ == "__main__":
    seterr(all="raise")
    grid = arange(1., 3250., 0.01)
    for molecule in ["CFC-11", "CFC-12", "CFC-113", "HCFC-123"]:
        mol = HitranCrossSection(molecule)
        for pressure in [10000., 100000.]:
            for temperature in [211.11, 250.23, 290.11]:
                xsec = mol.absorption_coefficient(temperature, pressure, grid)
                plt.plot(grid, xsec, label="{} K, {} Pa".format(temperature, pressure))
        plt.title(molecule)
        plt.legend()
