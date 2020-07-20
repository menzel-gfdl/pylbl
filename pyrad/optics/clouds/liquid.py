from netCDF4 import Dataset
from numpy import append, copy, insert, power, searchsorted, zeros

from .utils import CloudOptics
from ..utils import interp, Optics


class LiquidCloudOptics(CloudOptics):
    """Liquid water cloud optics parameterization.
       doi: 10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2

    Attributes:
        a1: Exctinction coefficient parameter.
        a2: Single-scatter albedo parameter.
        a3: Asymmetry factor parameter.
        band_limits: Lower/upper bounds of parameterization [cm-1].
        bands: Parameterization band centers [cm-1].
        b1: Exctinction coefficient parameter.
        b2: Single-scatter albedo parameter.
        b3: Asymmetry factor parameter.
        c1: Exctinction coefficient parameter.
        c2: Single-scatter albedo parameter.
        c3: Asymmetry factor parameter.
        max_radius: Maximum radius defined in parameterization [micron].
        min_radius: Minimum radius defined in parameterization [micron].
        radii: Radius bins [micron] for parameterization.
    """

    def __init__(self, path):
        with Dataset(path, "r") as dataset:
            radius = dataset.variables["radius"]
            self.min_radius, self.max_radius = radius.getncattr("valid_range")
            radius_bounds = dataset.variables[radius.getncattr("bounds")]
            self.radii = append(radius_bounds[:,0], radius_bounds[-1,-1])
            band = dataset.variables["band"]
            self.band_limits = band.getncattr("valid_range")
            self.bands = copy(band)
            for name in ("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3"):
                setattr(self, name, copy(dataset.variables[name]))

    def optics(self, water_concentration, equivalent_radius, grid):
        """Calculates cloud optics.

        Args:
            water_concentration: Water concentration [g m-3].
            equivalent_radius: Droplet equivalent radius [micron].
            grid: Spectral grid [cm-1].

        Returns:
            extinction_coefficient: Extinction coefficient [cm-1] (grid).
            single_scatter_albedo: Single-scatter albedo (grid).
            asymmetry_factor: Asymmetry factor (grid).
        """
        r = max(self.min_radius, equivalent_radius) if equivalent_radius < self.max_radius \
            else min(self.max_radius, equivalent_radius)
        i = searchsorted(self.radii, r) - 1

        n = self.bands.size + 1
        beta, omega, g = zeros(n + 1), zeros(n + 1), zeros(n + 1)
        cm_to_km = 1.e-5
        beta[1:n] = water_concentration*cm_to_km*(self.a1[i,:]*power(r, self.b1[i,:]) +
                    self.c1[i,:]) #Equation 13.
        omega[1:n] = 1. - (self.a2[i,:]*power(r, self.b2[i,:]) + self.c2[i,:]) #Equation 14.
        g[1:n] = self.a3[i,:]*power(r, self.b3[i,:]) + self.c3[i,:] #Equation 15.
        beta[0], omega[0], g[0] = beta[1], omega[1], g[1]
        beta[-1], omega[-1], g[-1] = beta[-2], omega[-2], g[-2]

        optics_ = Optics(grid)
        bands = append(insert(self.bands, 0, self.band_limits[0]), self.band_limits[-1])
        optics_.tau = interp(bands, beta, grid.points)
        optics_.omega = interp(bands, omega, grid.points)
        optics_.g = interp(bands, g, grid.points)
        return optics_
