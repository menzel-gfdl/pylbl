from pyrad.optics.gas import Gas
from pyrad.utils.grids import UniformGrid1D
import matplotlib.pyplot as plt

from pyrad.lbl.hitran import Voigt

grid = UniformGrid1D(600., 900., 0.01)
gas = Gas(formula="CO2", line_profile=Voigt())
k = gas.absorption_coefficient(temperature=299.7, pressure=101300.,
                               volume_mixing_ratio=.02595108, spectral_grid=grid.points)
plt.plot(grid.points, k)
plt.title("CO2 Absorption Spectrum")
plt.xlabel("wavenumber [cm-1]")
plt.ylabel("absorption coefficient [m2].")
plt.savefig("gas-optics.png")
