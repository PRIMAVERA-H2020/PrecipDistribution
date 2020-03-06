import iris
import add_laea
from iris.experimental.regrid_conservative import regrid_conservative_via_esmpy as rcve


def main():
 indir = '/gws/nopw/j04/primavera3/cache/sberthou/UKCPOBS/'
 outdir = '/gws/nopw/j04/primavera3/cache/sberthou/ukcpobs_regrid_CORDEX50/'
 target_grid =  iris.load_cube('/home/users/sberthou/PrecipExtremes/EUROCORDEX_grid.nc')
 cube_xy_guessbounds(target_grid)
 
 for year in range(1977, 1979):
   for month in range(1, 12):
      cube = iris.load_cube('{}/daily_rainfall_{}.nc'.format(indir, year))
      cube2 = rcve(cube, target_grid)
      iris.save(cube2, '{}/ukcpobs_regrid_CORDEX50_{}.nc'.format(outdir, year))


def cube_primary_xy_coord_names(cube):
    """Return the primary latitude and longitude coordinate standard names, or
    long names, from a cube.

    Arguments:
        cube (:class:`iris.cube.Cube`): An Iris cube

    Returns:
        The names of the primary latitude and longitude coordinates
    """
    check_cube_instance(cube)
    latc = cube.coords(axis='y')[0] if cube.coords(axis='y') else -1
    lonc = cube.coords(axis='x')[0] if cube.coords(axis='x') else -1

    if -1 in (latc, lonc):
        msg = "Error retrieving xy dimensions in cube: {!r}"
        raise ValueError(msg.format(cube))

    latitude = latc.standard_name if latc.standard_name else latc.long_name
    longitude = lonc.standard_name if lonc.standard_name else lonc.long_name
    return latitude, longitude

def check_cube_instance(cube):
    """Check an iris.Cube instance has been provided.

    Arguments:
        cube (:class:`iris.cube.Cube`): The cube to check

    Returns:
        `True` if the passed argument is a cube, `False` otherwise

    Raises:
        TypeError: If cube passed is not an Iris cube
    """
    if not isinstance(cube, iris.cube.Cube):
        msg = "Iris.Cube instance required, got {}"
        raise TypeError(msg.format(type(cube)))
    return True

def cube_xy_guessbounds(cube):
    """Guess latitude/longitude bounds of the cube and add them (**in place**)
    if not present.

    Arguments:
        cube (:class:`iris.cube.Cube`): An Iris cube

    Warning:
        This function modifies the passed `cube` in place, adding bounds to the
        latitude and longitude coordinates.
    """
    check_cube_instance(cube)
    for coord in cube_primary_xy_coord_names(cube):
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()

if __name__ == '__main__':
   main()
