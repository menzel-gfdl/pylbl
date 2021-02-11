from csv import reader
from ftplib import FTP
from os.path import join
from tarfile import TarFile
from tempfile import TemporaryDirectory

from numpy import asarray
from scipy.interpolate import interp1d


ADDRESS = "ftp2.gfdl.noaa.gov"
TARFILE = "/perm/GFDL_pubrelease/test_data/grtcode-data.tar.gz"
Pa_to_atm = 9.86923e-6
cm_to_m = 0.01


class GriddedField(object):
    """Container for one-dimensional data fields at the specified grid points.

    Attributes:
        data: Numpy array of data values.
        grid: Numpy array of grid points.
    """
    def __init__(self, data, grid):
        self.data = data
        self.grid = grid

    def regrid(self, grid, fill=0.):
        """Interpolates data to a new grid points.

        Args:
            grid: Numpy array of new grid points.

        Returns:
            Numpy array of interpolated data at the input grid points.
        """
        return interp1d(self.grid, self.data, bounds_error=False, fill_value=fill)(grid)


def download_gfdl_data(continuum, directory, names, parameters):
    """Downloads data from GFDL's FTP site."""
    with TemporaryDirectory() as tmp:
        archive = join(tmp, "archive.tar.gz")
        with open(archive, "wb") as tarball:
            with FTP(ADDRESS) as ftp:
                ftp.login()
                ftp.retrbinary("RETR {}".format(TARFILE), tarball.write)
        tarball = TarFile.open(name=archive, mode="r")
        for name, parameter in zip(names, parameters):
            path = join("grtcode-data", directory, name)
            tarball.extract(path, path=tmp)
            with open(join(tmp, path), newline="") as csvfile:
                next(csvfile)
                grid, data = [], []
                for row in reader(csvfile, delimiter=","):
                    grid.append(float(row[0]))
                    data.append(float(row[1]))
            setattr(continuum, parameter, GriddedField(asarray(data), asarray(grid)))
