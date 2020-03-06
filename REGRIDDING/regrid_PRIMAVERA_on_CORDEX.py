import european_masked_subregion as ems
import iris
import subprocess
import os
from cf_units import Unit

def main():
    outdir = '/gws/nopw/j04/primavera3/cache/sberthou/'
    country = 'prudence'
    frequency = 'd'
    n512 = False
    other = 'PRIMAVERA'
    time_con = iris.Constraint(time=lambda cell: 1950 <= cell.point.year < 2006)
    # CORDEX grid, obtained by executing "cdo griddes CORDEX_file.nc > EUROCORDEX_grid.txt"
    region_name = 'EUROCORDEX'
    new_grid = '{}_grid.txt'.format(region_name)

    runlist, _ = ems.get_runlist_region(frequency, n512, country, other=other)
    for model in runlist:
      #### 1st step: load all the data into one single cube and save it, so that CDO regrid is fast (calculates the weights only once, then applies it to the whole time series. ####
      #### at the very end, once you're happy that everything has been regridded, don't forget to delete the large files, as they are enormous! ####    
      large_cube = '{}/pr_{}'.format(outdir, model)
          
      if not os.path.exists('{}_with_grid.nc'.format(large_cube)):
         if not os.path.exists('{}.nc'.format(large_cube)):
            cubelist = iris.load(runlist[model], time_con, callback=callback_overwrite)
            iris.util.unify_time_units(cubelist)
            cube = cubelist.concatenate_cube()
            print('cubes loaded')
            iris.save(cube, large_cube+'.nc')
            print('large cube saved')
         elif 'EC-Earth' in model:
            cube = iris.load_cube('{}.nc'.format(large_cube))
            redefine_spatial_coords(cube)
            iris.save(cube, '{}_tmp.nc'.format(large_cube))
            cmd = 'mv {}_tmp.nc {}.nc'.format(large_cube, large_cube)
            shellcmd(cmd, 'mv ECarth failed')
      #### get the grid from the large netCDF file ####
         cmd = 'cdo griddes {}.nc > init_grid.txt'.format(large_cube)
         shellcmd(cmd, 'cdo griddes didn''t complete')
      #### set the grid in the file (not sure why you have to do this, I think you need it mostly for the grids which have 1D lat/lon, because cdo requires 2D lat,lons, calculated with cdo griddes ####
         cmd = 'cdo -setgrid,init_grid.txt {}.nc {}_with_grid.nc'.format(large_cube, large_cube)
         shellcmd(cmd, 'cdo setgrid didn''t complete')
      #### remapping itself ####
      if not os.path.exists('{}_regridded_on_{}.nc'.format(large_cube, region_name)):
        if not model == 'EC-Earth3P-HR':
          cmd = 'cdo remapcon,{} {}_with_grid.nc {}_regridded_on_{}.nc'.format(new_grid, large_cube, large_cube, region_name)
          shellcmd(cmd, 'cdo remapcon didn''t complete')
        else:
          cmd = 'cdo remapcon,{} {}.nc {}_regridded_on_{}.nc'.format(new_grid, large_cube, large_cube, region_name)
          shellcmd(cmd, 'cdo remapcon didn''t complete')

def shellcmd(cmd, msg):
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print('syst.cmd terminated by signal', retcode)
        elif retcode:
            print('syst.cmd returned in ', msg, '', retcode)
    except OSError as ex:
        print("Execution failed in " + msg + ": ", ex)


def callback_overwrite(cube, field, filename):
    coord2rm = ['forecast_reference_time', 'forecast_period', 'season_number',
                '3hours', 'hours']
    for co2rm in coord2rm:
        if co2rm in [coord.name() for coord in cube.coords()]:
            cube.remove_coord(co2rm)
    attributes_to_overwrite = ['date_created', 'log', 'converter', 'um_streamid',
                               'creation_date', 'history', 'iris_version', 'prod_date',
                               'CDI', 'CDO', 'ArchiveMetadata.0', 'CoreMetadata.0', 'creation_date', 'tracking_id' ]
    for att in attributes_to_overwrite:
        if cube.attributes.has_key(att):
            cube.attributes[att] = 'overwritten'
    attributes_to_del = ['radar.flags', 'log', 'iris_version', '_NCProperties', 'NCO']
    for att in attributes_to_del:
        if cube.attributes.has_key(att):
            del cube.attributes[att]

    if cube.coords('T'):  # for GPCP
        cube.coord('T').standard_name = 'time'


def redefine_spatial_coords(cube):
    """
    Redefines the latitude and longitude points for the EC-Earth3 model
    into single, rather than multi-dimensional, coordinates.
    """
    # procedure for handling EC-Earth latitude conversion
    cube.coord('cell index along second dimension').points = cube.coord(
      'latitude').points[:,0]
    cube.remove_coord('latitude') # remove AuxCoord 'latitude'
    cube.coord('cell index along second dimension') \
        .rename('latitude') # assign DimCoord 'latitude'
    cube.coord('latitude').units = Unit('degrees')
    cube.coord('latitude').long_name = 'latitude'
    cube.coord('latitude').var_name = 'lat'
    cube.coord('latitude').guess_bounds()
    # procedure for handling EC-earth longitude conversion
    cube.coord('cell index along first dimension').points = cube.coord(
       'longitude').points[0,:]
    cube.remove_coord('longitude') # remove AuxCoord 'longitude'
    cube.coord('cell index along first dimension') \
        .rename('longitude') # assign DimCoord 'longitude'
    cube.coord('longitude').units = Unit('degrees')
    cube.coord('longitude').long_name = 'longitude'
    cube.coord('longitude').var_name = 'lon'
    cube.coord('longitude').guess_bounds()


if __name__ == '__main__':
    main()
