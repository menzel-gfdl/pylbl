from numpy import append, insert, zeros
from scipy.interpolate import interp1d

from ..utils.grids import GridError


class Optics(object):
    def __init__(self, grid):
        self.grid = grid
        self.tau, self.omega, self.g = zeros(grid.size), zeros(grid.size), zeros(grid.size)

    def __add__(self, other):
        if self.grid != other.grid:
            raise GridError("grids do not match.")
        new = Optics(self.grid)
        new.tau[:] = self.tau[:] + other.tau[:]
        new.omega = (self.tau[:]*self.omega[:] + other.tau[:]*other.omega[:])/new.tau[:]
        new.g = (self.tau[:]*self.omega[:]*self.g[:] + other.tau[:]*other.omega[:]*other.g[:])/ \
                (new.tau[:]*new.omega[:])
        return new

    def __radd__(self, other):
        return self + other

def interp(x, y, newx):
    f = interp1d(x, y, bounds_error=False, fill_value=(y[0], y[-1]))
    return f(newx)
