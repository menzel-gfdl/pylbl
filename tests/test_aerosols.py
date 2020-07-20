from os.path import join

from pyrad.optics.aerosols import AerosolOptics
from pyrad.utils.grids import UniformGrid1D


grid = UniformGrid1D(1., 3000., .1)
concentration = 1. #
relative_humidity = 50.5 #[%].
mixture = 80.3 #[%].

sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
optics = sulfate.optics(concentration, grid, relative_humidity, mixture)
