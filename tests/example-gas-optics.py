from pyrad.optics.gas import Gas
from pyrad.utils.grids import UniformGrid1D
import matplotlib.pyplot as plt

spectral_grid = UniformGrid1D(1., 500., 0.1)
gas = Gas(formula="H2O")
k = gas.absorption_coefficient(temperature=300., pressure=1., partial_pressure=0.9,
                               spectral_grid=spectral_grid.points)

plt.plot(spectral_grid.points, k)
plt.title("H2O Absorption Spectrum")
plt.xlabel("wavenumber [cm-1]")
plt.ylabel("absorption coefficient [cm-1].")
plt.savefig("gas-optics.png")
