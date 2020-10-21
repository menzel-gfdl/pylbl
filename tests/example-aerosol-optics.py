from pyrad.optics.aerosols import AerosolOptics

sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
aerosol_optics = sulfate.optics(concentration, spectral_grid, relative_humidity, mixture)
