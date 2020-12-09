from netCDF4 import Dataset
from numpy import append, asarray, copy, insert, matmul, searchsorted, zeros

from .utils import CloudOptics
from ..utils import interp, Optics


class IceCloudOptics(CloudOptics):
    """Ice water cloud optics parameterization.
       doi: 10.1175/2009JCLI2844.1

    Attributes:
        a: a parameters from equations 4a/5a.
        b: b parameters from equations 4b/5b.
        bands: Parameterization band limits [cm-1].
        c: c parameters from equations 4c/5c.
        last_ir_band: Index of last infrared (longwave) band.
        radii: Radius bins [micron] for the parameterization.
    """

    def __init__(self, path):
        with Dataset(path, "r") as dataset:
            self.radii = copy(dataset.variables["radius_bnds"][...])
            band = dataset.variables["band_bnds"]
            self.bands = copy(band[...])
            self.last_ir_band = band.getncattr("last_IR_band")
            for name in ("a", "b", "c"):
                setattr(self, name, copy(dataset.variables[name]))

    def optics(self, ice_content, equivalent_radius, grid, mode="longwave"):
        """Calculates cloud optics.

        Args:
            ice_content: Ice content [g m-3].
            equivalent_radius: Particle equivalent radius [micron].
            grid: Spectral grid [cm-1].
            mode: Spectral region (longwave or shortwave) to consider.
        Returns:
            extinction_coefficient: Extinction coefficient [cm-1] (grid).
            single_scatter_albedo: Single-scatter albedo (grid).
            asymmetry_factor: Asymmetry factor (grid).
        """
        r = searchsorted(self.radii[:, 0], equivalent_radius) - 1
        i = self.last_ir_band
        d = asarray([equivalent_radius**i for i in range(self.a.shape[-1])])
        d_inv = 1./d[:]

        if mode.lower() == "longwave":
            n = i + 1
            bands = append(insert(0.5*(self.bands[:i, 0] + self.bands[:i, 1]), 0, self.bands[0, 0]),
                           self.bands[i-1, 1])
            a, b, c = self.a[:i, :], self.b[:i, :], self.c[r, :i, :]
        elif mode.lower() == "shortwave":
            n = self.bands.shape[0] - i + 1
            bands = append(insert(0.5*(self.bands[i:, 0] + self.bands[i:, 1]), 0, self.bands[i, 0]),
                           self.bands[-1, 1])
            a, b, c = self.a[i:, :], self.b[i:, :], self.c[r, i:, :]
        else:
            raise ValueError("mode must be either 'longwave' or 'shortwave'.")

        tau, omega, g = zeros(n+1), zeros(n+1), zeros(n+1)
        tau[1:n] = ice_content*(matmul(a, d_inv))
        if mode.lower() == "longwave":
            omega[1:n] = ice_content*(matmul(b, d_inv))
        else:
            omega[1:n] = 1. - matmul(b, d)
        g[1:n] = matmul(c, d)
        tau[0], omega[0], g[0] = tau[1], omega[1], g[1]
        tau[-1], omega[-1], g[-1] = tau[-2], omega[-2], g[-2]

        optics_ = Optics(grid)
        optics_.tau = interp(bands, tau, grid.points)
        optics_.omega = interp(bands, omega, grid.points)
        optics_.g = interp(bands, g, grid.points)
        return optics_
