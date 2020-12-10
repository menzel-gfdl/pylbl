from numpy import ceil, linspace


class GridError(BaseException):
    pass


class UniformGrid1D(object):
    """Uniform one-dimensional grid.

    Attributes:
        points: Numpy array of grid points.
    """

    def __init__(self, lower_bound, upper_bound, resolution):
        """Initializes a uniform one-dimensional grid.

        Args:
            lower_bound: Grid lower bound.
            upper_bound: Grid upper bound.
            resolution: Grid resolution.
        """
        size = int(ceil((upper_bound - lower_bound)/resolution) + 1)
        self.points = linspace(lower_bound, upper_bound, size, endpoint=True)

    @property
    def size(self):
        return self.points.size

    @property
    def start(self):
        return self.points[0]

    @property
    def stop(self):
        return self.points[-1]
