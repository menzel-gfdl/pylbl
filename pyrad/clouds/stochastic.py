from numpy import exp, nonzero, zeros
from numpy.ma import masked_where
from numpy.random import rand
from scipy.special import betainc, betaincinv


class TotalWaterPDF(object):
    """Tota water mixing ratio in dry air probability distribution function.

    Attributes:
        p: incomplete beta distribution shape parameter.
        q: incomplete beta distribution shape parameter.
    """

    def __init__(self, p, q):
        self.p = p
        self.q = q

    def specific_saturation_humidity(self, cloud_fraction):
        """Calculates normalized saturation specific humidity from equation A1 from
           doi: 10.1175/MWR3257.1

        Args:
            cloud_fraction: saturated volume fraction.

        Returns:
            normalized saturation specific humidity [kg/kg].
        """
        return betaincinv(self.p, self.q, 1. - cloud_fraction)

    def width(self, cloud_fraction, lwc, iwc):
        """Calculates the width of the total water probability distribution function (b - a)
           from equation A2 from doi: 10.1175/MWR3257.1, ignoring the parameter alpha.

        Args:
            cloud_fraction: saturated volume fraction.
            lwc: cloud liquid water condensate mixing ratio in dry air [kg/kg].
            iwc: cloud ice water condensate mixing ratio in dry air [kg/kg].

        Returns:
            distribution width (b - a).
        """
        qs = self.specific_saturation_humidity(cloud_fraction)
        return (lwc + iwc)/((self.p/(self.p + self.q))*(1. - betainc(self.p + 1, self.q, qs)) \
               - qs*cloud_fraction)

    def sample_condensate(self, cloud_fraction, lwc, iwc, overlap):
        """Calculates liquid and ice condensate amounts in each layer, based on the
           description from the appendix of doi: 10.1175/MWR3257.1.

        Args:
            cloud_fraction: saturated volume fraction.
            lwc: cloud liquid water condensate mixing ratio in dry air [kg/kg].
            iwc: cloud ice water condensate mixing ratio in dry air [kg/kg].

        Returns:
            Masked arrays of liquid and ice condensate amounts in each layer.
        """
        x = cloudiness(cloud_fraction, overlap)
        qa = masked_where(x.mask, cloud_fraction)
        qs = self.specific_saturation_humidity(qa)
        liquid, ice = masked_where(x.mask, lwc), masked_where(x.mask, iwc)
        width = self.width(qa, liquid, ice)
        total_condensate = width*(betaincinv(self.p, self.q, x) - qs)
        liquid_fraction = liquid/(liquid + ice)
        return total_condensate*liquid_fraction, total_condensate*(1. - liquid_fraction)


def overlap_parameter(altitude, scale_length):
    """Calculates the overlap parameter defined in equation 2 from doi: 10.1029/2004JD005100.

    Args:
        altitude: layer altitudes, ordered from top of atmosphere to the surface.
        scale_length: scale length parameter.

    Returns:
        overlap parameter between adjacent layers.
    """
    return exp((altitude[1:] - altitude[:-1])/scale_length)


def cloudiness(cloud_fraction, overlap_parameter):
    """Calculates boolean array for cloudiness from equation 1 from doi: 10.1256/qj.03.99.

    Args:
        cloud_fraction: layer saturation volume fraction.
        overlap_paramter: overlap parameter between adjacent layers.

    Returns:
        Masked rank array of random numbers, where mask=False values indicate the
        presence of a cloud.
    """
    x, r = rand(cloud_fraction.size), rand(cloud_fraction.size - 1)
    for i in nonzero(r <= overlap_parameter)[0]:
        x[i+1] = x[i]
    return masked_where(x <= 1. - cloud_fraction, x)
