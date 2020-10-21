from os.path import join
from pyrad.optics.clouds.ice import IceCloudOptics
from pyrad.optics.clouds.liquid import LiquidCloudOptics

ice_clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
ice_particle_size = 10. #microns.
ice_cloud_optics = ice_clouds.optics(iwc, ice_particle_size, spectral_grid,
                                     mode="longwave")
liquid_clouds = LiquidCloudOptics(join("pyrad_data", "clouds", "hu_stamnes.nc"))
liquid_droplet_radius = 10. #microns.
liquid_cloud_optics = liquid_clouds.optics(lwc, liquid_droplet_radius, spectral_grid)
