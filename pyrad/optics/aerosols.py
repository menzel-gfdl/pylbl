from netCDF4 import Dataset
from numpy import asarray, copy, reshape, searchsorted, zeros

from .utils import interp, Optics


class AerosolOptics(object):
    def __init__(self, path):
        with Dataset(path, "r") as dataset:
            self.bands = copy(dataset.variables["wavenumber"][:])
            try:
                self.humidity, self.humidity_map = self._make_map(dataset, "relative_humidity")
            except KeyError:
                self.humidity, self.humidity_map = None, None
            try:
                self.mixture, self.mixture_map = self._make_map(dataset, "mixture")
            except KeyError:
                self.mixture, self.mixture_map = None, None
            beta = dataset.variables["extinction_coefficient"]
            self.name = beta.getncattr("species")
            self.extinction_coefficient = copy(beta[...])
            self.single_scatter_albedo = copy(dataset.variables["single_scatter_albedo"][...])
            self.asymmetry_factor = copy(dataset.variables["asymmetry_factor"][...])

    @staticmethod
    def _make_map(dataset, name):
        x = copy(dataset.variables[name][:])
        x_map = asarray([searchsorted(x, i) for i in range(101)])
        return x, x_map

    def optics(self, concentration, grid, humidity=0, mixture=0):
        """

        Args:
            concentration: Concentration [kg m-2].
            humidity: Relative humidity percentage.
            mixture: Mixture percentage.

        Returns:
            An Optics object.
        """
        g_per_kg = 1.e3 #![g kg-1]
        i, x = (0, 1) if self.humidity_map is None else \
               (self.humidity_map[int(humidity)], self.humidity.size)
        j, y = (0, 1) if self.mixture_map is None else \
               (self.mixture_map[int(mixture)], self.mixture.size)
        shape = (y, x, self.extinction_coefficient.size//(x*y))

        tau = reshape(self.extinction_coefficient, shape)[j,i,:]*concentration*g_per_kg
        omega = reshape(self.single_scatter_albedo, shape)[j,i,:]
        g = reshape(self.asymmetry_factor, shape)[j,i,:]

        optics_ = Optics(grid)
        optics_.tau = interp(self.bands, tau, grid.points)
        optics_.omega = interp(self.bands, omega, grid.points)
        optics_.g = interp(self.bands, g, grid.points)
        return optics_
