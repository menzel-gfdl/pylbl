pyrad
=====

.. image:: https://travis-ci.com/menzel-gfdl/pylbl.svg?branch=reorganize-take-2
   :target: https://travis-ci.com/menzel-gfdl/pylbl
   :alt: Build Status

.. image:: https://readthedocs.org/projects/pylbl/badge/?version=latest
   :target: https://pylbl.readthedocs.io/en/latest/
   :alt: Documentation Status

Pyrad is a simple (one-dimensional), pure python, all-sky atmospheric radiation package.


Gas optics are calculated in a line-by-line fashion, using data from the HITRAN database:

.. code-block:: python

   from pyrad.lbl.hitran.database import Database
   from pyrad.lbl.line_profiles import Voigt
   from pyrad.optics.gas import Gas
   from pyrad.utils.grids import UniformGrid1D

   spectral_grid = UniformGrid1D(1., 3000., 0.1)
   hitran = Database()
   gas = Gas(formula, hitran)
   k = gas.absorption_coefficient(temperature, pressure, partial_pressure,
                                  spectral_grid.points, Voigt())

Clouds are generated in a stochastic fashion (typically found in GCMs):

.. code-block:: python

   from pyrad.optics.clouds.stochastic import overlap_parameter, TotalWaterPDF

   overlap = overlap_parameter(altitude, scale_length)
   water = TotalWaterPDF()
   lwc, iwc = water.sample_condensate(cloud_fraction, liquid_water_content,
                                      ice_water_content, overlap)

and their optics are calculated using standard look-up table parameterizations:

.. code-block:: python

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

Aerosol optics are also calculated using a GCM parameterization:

.. code-block:: python

   from pyrad.optics.aerosols import AerosolOptics

   sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
   aerosol_optics = sulfate.optics(concentration, spectral_grid, relative_humidity, mixture)


