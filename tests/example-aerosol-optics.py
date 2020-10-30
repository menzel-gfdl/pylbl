from os.path import join
from pyrad.optics.aerosols import AerosolOptics
from pyrad.utils.grids import UniformGrid1D

sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
aerosol_optics = sulfate.optics(concentration=0.5, grid=UniformGrid1D(1., 500., 0.1),
                                relative_humidity=50, mixture=50)
