from numpy import append, insert
from scipy.interpolate import interp1d


class CloudOptics(object):
    def optics(water, equivalent_radius, grid):
        raise NotImplementedError("this must be overridden.")


def interp(x, xlims, y, newx):
    z = append(insert(x, 0, xlims[0]), xlims[-1])
    f = interp1d(z, y, bounds_error=False, fill_value=(y[0], y[-1]))
    return f(newx)
