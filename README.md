pyrad
=====

[![Documentation status](https://readthedocs.org/projects/pylbl/badge/?version=latest)](https://pylbl.readthedocs.io/en/latest/)

Pyrad is a simple (one-dimensional), pure python, all-sky atmospheric radiation package.

Gas absorption coefficients can be calculated by:

```python
from pyrad.lbl.hitran import Voigt
from pyrad.optics.gas import Gas
from pyrad.utils.grids import UniformGrid1D

gas = Gas(formula="CO2", line_profile=Voigt())
k = gas.absorption_coefficient(temperature=299.7, pressure=101300., volume_mixing_ratio=.02595108,
                               spectral_grid=UniformGrid1D(600., 900., 0.01).points)
```

By default, the above code will download the necessary molecular line and total partition
function data from the web.  This can take a significant amount of time, especially if the
Gas objects are created often.  To retify this, I recommend creating local SQLite databases,
then re-using when creating Gas objects:

```python
from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.optics.gas import Gas

for x in ["H2O", "CO2", "O3"]:
    Hitran(x, Voigt()).create_database("hitran.sqlite")
    TotalPartitionFunction(x).create_database("tips-2017.sqlite")
gas = Gas("H2O", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
gas = Gas("CO2", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
gas = Gas("O3", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
```

Clouds are generated in a stochastic fashion (typically found in GCMs):

```python
from numpy import array, zeros
from pyrad.optics.clouds.stochastic import overlap_parameter, TotalWaterPDF

altitude = array([1., 2., 3.])
cloud_fraction = array([.9, .8, .75])
liquid_water_content = array([.3, .3, .3])
ice_water_content = array([.2, .2, .2])
num_subcolumns = 10
lwc, iwc = zeros((num_subcolumns, altitude.size)), zeros((num_subcolumns, altitude.size))

overlap = overlap_parameter(altitude, scale_length=2.)
for i in range(num_subcolumns):
    # Repeat sampling in each column.
    lwc[i, :], iwc[i, :] = TotalWaterPDF().sample_condensate(cloud_fraction, liquid_water_content,
                                                             ice_water_content, overlap=overlap)
```

and their optics are calculated using standard look-up table parameterizations:

```python
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
```
                                           
Aerosol optics are also calculated using a GCM parameterization:

```python
from os.path import join
from pyrad.optics.aerosols import AerosolOptics
from pyrad.utils.grids import UniformGrid1D

sulfate = AerosolOptics(join("pyrad_data", "aerosols", "sulfate_optics.nc"))
aerosol_optics = sulfate.optics(concentration=0.5, grid=UniformGrid1D(1., 500., 0.1),
                                relative_humidity=50, mixture=50)
```
