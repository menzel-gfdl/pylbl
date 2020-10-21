from pyrad.optics.clouds.stochastic import overlap_parameter, TotalWaterPDF
import matplotlib.pyplot as plt
from numpy import array, zeros

altitude = array([1., 2., 3.,])
cloud_fraction = array([.9, .8, .75])
liquid_water_content = array([.3, .3, .3])
ice_water_content = array([.2, .2, .2])
num_subcolumns = 10
lwc, iwc = zeros((num_subcolumns, altitude.size)), zeros((num_subcolumns, altitude.size))

overlap = overlap_parameter(altitude, scale_length=2.)
for i in range(num_subcolumns):
    #Repeat sampling in each column.
    lwc[i,:], iwc[i,:] = TotalWaterPDF().sample_condensate(cloud_fraction, liquid_water_content,
                                                           ice_water_content, overlap=overlap)

for condensate, name in zip([lwc, iwc], ["Liquid", "Ice"]):
    plt.clf()
    plt.imshow(condensate.transpose(), cmap="Greys", interpolation="none")
    plt.title("{} water content".format(name))
    plt.xlabel("sub-column")
    plt.xticks([], [])
    plt.ylabel("layer")
    plt.yticks([], [])
    colorbar = plt.colorbar()
    colorbar.set_label("g m-2")
    plt.savefig("{}-stochastic-clouds.png".format(name).lower())
