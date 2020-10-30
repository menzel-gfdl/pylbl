from logging import basicConfig, INFO, info
from os.path import join

from numpy import asarray, seterr

from pyrad.lbl.hitran import Voigt
from pyrad.optics.aerosols import AerosolOptics
from pyrad.optics.gas import Gas
from pyrad.optics.clouds.ice import IceCloudOptics
from pyrad.optics.clouds.liquid import LiquidCloudOptics
from pyrad.optics.clouds.stochastic import overlap_parameter, TotalWaterPDF
from pyrad.utils.grids import UniformGrid1D


if __name__ == "__main__":
    basicConfig(format="%(asctime)-15s - %(pathname)s(%(lineno)d):\n\t%(message)s",
                level=INFO)
    seterr(divide="raise", over="raise", invalid="raise")
    info("Starting integration test.")

    #Configure layer.
    mb_to_atm = 0.000986923
    spectral_grid = UniformGrid1D(1., 3000., 0.01)
    temperature = 288.99 #K
    pressure = 983.88*mb_to_atm #atm
    abundance = {"H2O" : 0.006637074, 
                 "CO2" : 0.0003599889,
                 "O3" : 6.859128e-08,
                 "N2O" : 3.199949e-07,
                 "CH4" : 1.700002e-06}

    #Gas optics (line-by-line).
    for formula, concentration in abundance.items():
        partial_pressure = pressure*concentration #atm

        #Calculate absorption coefficient.
        info("Starting absorption coefficient calculation.")
        gas = Gas(formula, line_profile=Voigt())
        k = gas.absorption_coefficient(temperature, pressure, partial_pressure,
                                       spectral_grid.points)
        info("Finished absorption coefficient calculation.")

    #Calculate cloud optics.
    cloud_fraction = asarray([0.99])
    ice_water_content =  asarray([50.]) #g m-3
    liquid_water_content = asarray([50.]) #g m-3
    altitude = asarray([5., 10.])
    scale_length = 2.
    overlap = overlap_parameter(altitude, scale_length)
    water = TotalWaterPDF()
    lwc, iwc = water.sample_condensate(cloud_fraction, liquid_water_content,
                                       ice_water_content, overlap)
    liquid_clouds = LiquidCloudOptics(join("tests", "pyrad_data", "clouds", "hu_stamnes.nc"))
    liquid_droplet_radius = 10. #micron
    liquid_cloud_optics = liquid_clouds.optics(lwc, liquid_droplet_radius, spectral_grid)
    ice_clouds = IceCloudOptics(join("tests", "pyrad_data", "clouds", "chou_suarez.nc"))
    ice_particle_size = 10. #micron
    ice_cloud_optics = ice_clouds.optics(iwc, ice_particle_size, spectral_grid,
                                         mode="longwave")

    #Calculate aerosol optics.
    concentration = 1. #
    relative_humidity = 50.5 #[%].
    mixture = 80.3 #[%].
    sulfate = AerosolOptics(join("tests", "pyrad_data", "aerosols", "sulfate_optics.nc"))
    aerosol_optics = sulfate.optics(concentration, spectral_grid, relative_humidity, mixture)

    info("Finishing integration test.")
