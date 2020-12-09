from os.path import join
from pyrad.optics.clouds.ice import IceCloudOptics
from pyrad.optics.clouds.liquid import LiquidCloudOptics
from pyrad.utils.grids import UniformGrid1D

ice_clouds = IceCloudOptics(join("pyrad_data", "clouds", "chou_suarez.nc"))
ice_cloud_optics = ice_clouds.optics(iwc=0.2, ice_particle_size=10,
                                     spectral_grid=UniformGrid1D(600., 900., .01),
                                     mode="longwave")

liquid_clouds = LiquidCloudOptics(join("pyrad_data", "clouds", "hu_stamnes.nc"))
liquid_cloud_optics = liquid_clouds.optics(lwc=0.2, liquid_droplet_radius=10.,
                                           spectral_grid=UniformGrid1D(500., 800., .1))
