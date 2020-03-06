import iris
from iris.experimental.regrid_conservative import regrid_conservative_via_esmpy as rcve
import iris.coord_categorisation as icat
import cf_units
import iris.analysis as iana

def main():
  indir = '/gws/nopw/j04/primavera3/cache/sberthou/FRANCE/'
  outdir = '/gws/nopw/j04/primavera3/cache/sberthou/FRANCE_regrid_CORDEX50/'
  target_grid =  iris.load_cube('/home/users/sberthou/PrecipExtremes/EUROCORDEX_grid.nc')
  cube_xy_guessbounds(target_grid)
  for year in range(1999, 2019):
      time_con = iris.Constraint(time=lambda cell: cell.point.year == year )
      #time_con1 = iris.Constraint(time=lambda cell: cell.point.month == 10 ) 
      yearm1 = year - 1
      cubelist = iris.cube.CubeList([])
      for var in ['PRCP', 'SNOW']:
         filey = '{}/Forc{}_france_SAFRAN_8Km_1hour_{}*_V2.nc'.format(indir, var, year)
         fileym1 = '{}/Forc{}_france_SAFRAN_8Km_1hour_{}*_V2.nc'.format(indir, var, yearm1)
         cubes = iris.load([filey, fileym1], 'precipitation_flux' & time_con, callback=callback_overwrite)
         iris.util.unify_time_units(cubes)
         hcube = cubes.concatenate_cube()
         hcube.long_name = None
         cubelist.append(hcube)
      hourly_cube = cubelist[0] + cubelist[1]
      cube = hourly2daily(hourly_cube)
      add_grid_to_safran(cube)
      cube2 = rcve(cube, target_grid)
      iris.save(cube2, '{}/ForcPRCPT_france_SAFRAN_CORDEX50_daily_{}.nc'.format(outdir, year))


def callback_overwrite(cube, field, filename):

    attributes_to_overwrite = ['actual_range', 'original_name']
    for att in attributes_to_overwrite:
        if cube.attributes.has_key(att):
            cube.attributes[att] = 'overwritten'
    
    if cube.coords('time'):
        cube.coord('time').attributes = None


def add_grid_to_safran(cube):
    grid_info = iris.load_cube('safran_grid.nc')
    lat_coord = grid_info.coord('projection_y_coordinate')
    lon_coord = grid_info.coord('projection_x_coordinate')
    #coord_sys = iris.coord_systems.LambertConformal(central_lat=grid_info.attributes['latitude_of_projection_origin'], 
    #                                            central_lon=grid_info.attributes['longitude_of_central_meridian'], 
    #                                            false_northing=grid_info.attributes['false_northing'], 
    #                                            false_easting=grid_info.attributes['false_easting'],
    #                                            grid_mapping_name = 'lambert_conformal_conic')
    # grid_tab_lat=range(2102500, 3202500, 5000) # missing info
    # grid_tab_lon = range(3502500, 4702500, 5000) # missing info

    #lat_coord = iris.coords.DimCoord(grid_tab_lat, standard_name = 'grid_latitude', 
    #                                            units=iunit.Unit('m'), 
    #                                            coord_system = coord_sys)
    #lon_coord = iris.coords.DimCoord(grid_tab_lon, standard_name = 'grid_longitude', 
    #                                            units=iunit.Unit('m'), 
    #                                            coord_system = coord_sys)
    cube.add_dim_coord(lat_coord, 1)
    cube.add_dim_coord(lon_coord, 2)
    cube.coord('projection_x_coordinate').units='m'
    cube.coord('projection_y_coordinate').units='m'


def hourly2daily(hourly_cube):
    '''
    Aggregates the hourly precipitation cube into daily_cube
    :param iris.cube: hourly cube to be aggregated in daily
    '''
    # raise exception is not only one year in a cube?
    icat.add_day_of_year(hourly_cube,  'time')
    hourly_cube.coord('time').bounds = None
    #print(hourly_cube.coord('day_of_year'))
    print(hourly_cube.coord('day_of_year').points[-80:])
    print(hourly_cube.coord('day_of_year').points[:80])
    #print(hourly_cube.coord('day_of_year').bounds)
    if hourly_cube.units in ['kg/m2/s', 'kg m-2 s-1', 'mm day-1', 'mm/day']:
        daily_data = hourly_cube.aggregated_by('day_of_year', iana.MEAN)
    elif hourly_cube.units in ['mm', 'mm/h', 'mm h-1', 'mm/3hr', 'mm hour-1']:
        print('summing hourly')
        daily_data = hourly_cube.aggregated_by('day_of_year', iana.SUM)
        daily_data.units = cf_units.as_unit('mm day-1')
    else:
        raise ValueError('wrong units: {}'.format(hourly_cube.units))
    # daily_data.name = 'precipitation_rate'
    return daily_data


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
