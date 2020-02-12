import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
# from matplotlib import style
import cartopy.crs as ccrs
from shapely.geometry.polygon import Polygon
from shapely.geometry.point import Point
import numpy as np
import iris
import os
import cPickle as pickle
import iris.analysis.cartography as iac
import json
import glob
from iris.cube import Cube
import iris.analysis as iana
import iris.coord_categorisation as icat
import iris.quickplot as qplt
from scipy import stats
import cf_units
import netCDF4
from cf_units import Unit
import iris.plot as iplt
import re
import cartopy.feature as cfeature
from matplotlib.colors import from_levels_and_colors
# from statsmodels.sandbox.stats import multicomp
from iris.analysis import maths as iam
import random
from mpl_toolkits.axes_grid1.inset_locator import inset_axes ###update11/10
from iris.coord_systems import LambertConformal

iris.FUTURE.cell_datetime_objects = True
font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}

matplotlib.rc('font', **font)
# plt.rcParams['text.usetex'] = True
# plt.rcParams['text.latex.unicode'] = True
RADIUS_EARTH = 6371229
SURF_AREA_EARTH = 4.0 * np.pi * RADIUS_EARTH ** 2


##### The dictionary hereunder is not used if regions are defined with a mask like process_data_europe.py ############
##### only used with process_data.py #####
srex_reg = {  ############# longitude must be between -180 and 180 ###############
    # 'WAF':  Polygon(((6.00, 5.00), (-13.0, 5.00), (-13.0, 17.5), (6.0, 17.5))),
    'GOG': Polygon(((5.00, 3.00), (-8.0, 3.00), (-8.0, 9.0), (5.0, 9.0))),
    'GSL': Polygon(((-8.0, 3.), (-18.0, 3.0), (-18.0, 12.0), (-8.0, 12.0))),
    'CSA': Polygon(((5.00, 9.0), (-8.00, 9.0), (-8.0, 15.0), (5.00, 15.0))),
    'NSA': Polygon(((5.00, 15.0), (-11.0, 15.0), (-11.0, 20.0), (5.00, 20.0))),
    'SEN': Polygon(((-11.00, 12.0), (-18.0, 12.0), (-18.0, 17.0), (-11.0, 17.0))),

    'NIN': Polygon(((-120.00, -7.5), (-178.0, -7.5), (-178.0, 7.5), (-120.0, 7.5))),
    'NWP': Polygon(((160.00, 5), (109.0, 10.0), (112.0, 20.0), (140.0, 35.0), (160.00, 5),)),
    # 'ALL': Polygon(((-18.0, 3.0), (-18.0, 20.0), (9.0, 20.0), (9.0, 3.0))),
                #'SAH':  Polygon(((25.0, 10.0), (-15.0, 10.0), (-15.0, 20.0), (25.00, 20.0))),
                'EAF':  Polygon(((34.5, 9.0), (29.0, 9.0), (29.0, 14.0), (34.5, 14.0))),
               'ALA':  Polygon(((-105.000, 60.00), (-168.022, 60.00), (-168.022, 72.55), (-105.000, 72.55))),
               'AMZ':  Polygon(((-66.377, -20.000), (-79.729, -1.239), (-68.800, 11.43), (-50.000, 11.43), (-50.000, -20.000))),
                'ANT':  Polygon(((-180.00, -90.01), (180.00, -90.01), (180.00, -60.00), (-180.00, -60.00))),
                'CAM':  Polygon(((-68.800, 11.43), (-79.729, -1.239), (-118.323, 28.56), (-90.315, 28.56))),
                'CAS':  Polygon(((60.00, 30.00), (60.00, 50.00), (75.00, 50.00), (75.00, 30.00))),
                'NEU':  Polygon(((-10.000, 48.00), (-10.000, 75.00), (40.00, 75.00), (40.00, 61.32))),
                'CEU':  Polygon(((-10.000, 45.00), (-10.000, 48.00), (40.00, 61.32), (40.00, 45.00))),
               'EUR':  Polygon(((-10.000, 45.00), (-10.000, 75.00), (40.00, 75.00), (40.00, 45.00))),
               'CGI':  Polygon(((-10.000, 50.00), (-105.000, 50.00), (-105.000, 85.00), (-10.000, 85.00))),
               'CNA':  Polygon(((-85.000, 50.00), (-85.000, 28.56), (-105.000, 28.56), (-105.000, 50.00))),
                'EAF':  Polygon(((25.00, -11.365), (25.00, 15.00), (51.99, 15.00), (51.99, -11.365))),
               'EAS':  Polygon(((100.00, 20.00), (100.00, 50.00), (145.00, 50.00), (145.00, 20.00))),
                'ENA':  Polygon(((-60.000, 25.00), (-85.000, 25.00), (-85.000, 50.00), (-60.000, 50.00))),
                'MED':  Polygon(((-10.000, 30.00), (-10.000, 45.00), (40.00, 45.00), (40.00, 30.00))),
               'NAS':  Polygon(((40.00, 50.00), (40.00, 70.00), (179.00, 70.00), (179.00, 50.00))),
                'NAU':  Polygon(((110.00, -30.000), (110.00, -10.000), (155.00, -10.000), (155.00, -30.000))),
                'NEB':  Polygon(((-34.000, -20.000), (-50.000, -20.000), (-50.000, 0.00), (-34.000, 0.00))),
                'NEU':  Polygon(((-10.000, 48.00), (-10.000, 75.00), (40.00, 75.00), (40.00, 61.32))),
                'SAF':  Polygon(((-10.000, -35.000), (-10.000, -11.365), (51.99, -11.365), (51.99, -35.000))),
                #'SAH':  Polygon(((-20.000, 15.00), (-20.000, 30.00), (40.00, 30.00), (40.00, 15.00))),
                'SAS':  Polygon(((60.00, 5.00), (60.00, 30.00), (100.00, 30.00), (100.00, 20.00), (95.00, 20.00), (95.00, 5.00))),
                'SAU':  Polygon(((110.00, -50.000), (110.00, -30.000), (180.00, -30.000), (180.00, -50.000))),
               'SEA':  Polygon(((95.00, -10.000), (95.00, 20.00), (155.00, 20.00), (155.00, -10.000))),
                'SSA':  Polygon(((-39.376, -20.000), (-39.376, -56.704), (-67.348, -56.704), (-72.141, -50.000), (-66.377, -20.000))),
                'TIB':  Polygon(((75.00, 30.00), (75.00, 50.00), (100.00, 50.00), (100.00, 30.00))),
               'WAF':  Polygon(((-20.000, -11.365), (-20.000, 15.00), (25.00, 15.00), (25.00, -11.365))),
                'WAS':  Polygon(((40.00, 15.00), (40.00, 50.00), (60.00, 50.00), (60.00, 15.00))),
               'WNA':  Polygon(((-105.000, 28.56), (-130.000, 28.56), (-130.000, 60.00), (-105.000, 60.00))),
                'WSA':  Polygon(((-79.729, -1.239), (-66.377, -20.000), (-72.141, -50.000), (-67.348, -56.704), (-82.022, -56.704), (-82.022, 0.53))),
##################### chap.10 AR6 regions ############################
               'SAH':  Polygon(((40.0, 10.0), (-15.0, 10.0), (-15.0, 20.0), (40.00, 20.0))), # JJAS
               'SESA':  Polygon(((-45.0, -40.0), (-65.0, -40.0), (-65.0, -25.0), (-45.00, -25.0))), #DJF
               'CAR':  Polygon(((-60.0, 11.0), (-85.0, 11.0), (-85.0, 25.0), (-60.00, 25.0))), # JJA
               'SWCA':  Polygon(((-94.0, 25.0), (-124.0, 25.0), (-124.0, 40.0), (-94.0, 40.0))), #std seasons?
               'SAU':  Polygon(((110.0, -40.0), (110.0, -25.0), (150.0, -25.0), (150.0, -40.000))),
               'EAS':  Polygon(((100.0, 20.0), (100.0, 50.0), (145.0, 50.00), (145.0, 20.00))), #JJA # need to check bounds, as diff btw 2 regions in chap10
			   'CEEU': Polygon(((140.0, 40.0), (60.0, 40.0), (60.0, 65.0), (140.00, 65.0))), # DJF (but temperature)
              # 'WEU': Polygon(((30.0, 35.0), (-15.0, 35.0), (-15.0, 70.0), (30.00, 70.0))) #probably don't need western Europe as already done in a lot of sub-regions
              # HIMA: to be defined + season?
}


class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return np.asscalar(obj)
        return json.JSONEncoder.default(self, obj)


def numpy_walk(node):
    assert isinstance(node, dict), 'This function only works on dictionaries'
    for key, item in node.items():
        if isinstance(item, list):
            node[key] = np.array(item)
        else:
            node[key] = item  # 0
        # numpy_walk(node)

########## this is useful if you have troubles reading files (concatenation error) ####
### you can add anything to remove attributes and so on ####
def callback_overwrite(cube, field, filename):
    coord2rm = ['forecast_reference_time', 'forecast_period', 'season_number',
                '3hours', 'hours']
    for co2rm in coord2rm:
        if co2rm in [coord.name() for coord in cube.coords()]:
            cube.remove_coord(co2rm)
    attributes_to_overwrite = ['date_created', 'log', 'converter', 'um_streamid',
                               'creation_date', 'history', 'iris_version', 'prod_date', 
                               'CDI', 'CDO', 'NCO', 'ArchiveMetadata.0', 'CoreMetadata.0', 
                               'tracking_id', 'starttime', 'endtime', 'run_startdate', 'actual_range', 'experiment' ]
    for att in attributes_to_overwrite:
        if cube.attributes.has_key(att):
            if att == 'history':
                if 'cdo remapcon' in cube.attributes[att]:
                    if 'CORDEX' in cube.attributes[att]:
                        cube.attributes['regridded_to'] = cube.attributes['history'].split(',')[1][:10]
            cube.attributes[att] = 'overwritten'
    attributes_to_del = ['radar.flags', 'log', 'iris_version', '_NCProperties', 'NCO', 'atmosphere_time_step_in_seconds', 'model_simulation_start', 'model_simulation_end', 'model_is_restarted']
    for att in attributes_to_del:
        if cube.attributes.has_key(att):
            del cube.attributes[att]

    if cube.coords('T'):  # for GPCP
        cube.coord('T').standard_name = 'time'


class PrecipExtremes(object):
    '''
    This class consists of an object composed of:
	- a region dictionary (name:Polygon)
        - a config dictionary, with infiles (), lsmask (mask file)
        - original_data: will be populated with the precipitation cube when load_mask_ppn is called
        - resolutions (str): resolution to be processed 
        - hists (dict of arrays): will be populated when calling gendata if hist is 1D (histogram data for every region at every resolution for every model run)
        - hist2d (dict of arrays): will be populated when calling gendata if hist is 2D (wet/dry spells)
        - bins (np.array): intensity bins
        - duration_bins (np.array): duration bins (only for 2D histograms)
        - totnum (dict of floats): will be populated with the number of spells when gendata is called for 2D histograms
        - start_year: start year in the processed file - this may be changed by changing the time_con in process_data.py
        - end_year: end year in the processed file - this may be changed by changing the time_con in process_data.py
        - nb_spatial_points: number of land points in the region of interest
        - centiles calculated on the histogram
        - yearly time series is just total rainfall per year, not used
        - drought_bins only used if wet/dry spells used, defined in process_data.py
    '''

    def __init__(self, srex_reg, config, resolution, bins=np.arange(200) / 86400., duration_bins=None,
                 drought_bins=None):
        self.regions = srex_reg
        self.config = config
        self.original_data = {}
        self.resolutions = resolution
        self.hists = {}
        self.hist2d = {}
        self.bins = bins
        self.duration_bins = duration_bins
        self.totnum = {}
        self.start_year = {}
        self.end_year = {}
        self.nb_spatial_points = {}
        self.centiles = {}
        self.yearly_time_series = {}
        self.drought_bins = drought_bins

    def load_mask_ppn(self, res, infile=None, maskfile=None, rain_and_snow=None, time_con=None, season=None,
                      season_list=('djf', 'mam', 'jja', 'son')):

        '''
	    Loads the mask and the precipitation file (ppn) for the data: the mask will be used to calculate the weights for the cube

        :param str res: String describing what the resolution of the data is 
        :param str infile: if infile present, loads it (them), otherwise, loads them from self.config['infiles']
        :param str maskfile: if a mask is provided from an outside file, loads it, otherwise, loads the mask from self.config['lsmasks'][res]
        :param bool rain_and_snow: if True, loads the CPM files and adds rain and snow together
        :param iris.Constraint time_con: restricts the CPM file loading to the time-constraint
        '''

        if rain_and_snow:
            self.original_data['ppn'] = loadpp_regridded2eobs(self.config['infiles'], time_con)
        else:
            if not infile:
                infile = self.config['infiles']  # [res]
            cubes = iris.load(infile, time_con, callback=callback_overwrite)
            season_cubes = iris.cube.CubeList([])
            # season = 'jas'
            # season_list = ('jfm', 'amj', 'jas', 'ond')
            if season:
                print 'extracting season: ', season
                for cube in cubes:
                    if not cube.coords('season'):
                        icat.add_season(cube, 'time', 'season', season_list)
                    cube.coord('time').bounds = None
                    season_con = iris.Constraint(season=lambda s: s == season)
                    season_cube = cube.extract(season_con)
                    if season_cube:
                        season_cubes.append(season_cube)
                cubes = season_cubes
            iris.util.unify_time_units(cubes)
            #### lots of conditions for specific datasets here ####
            if 'u-ad251' in infile[0]:
                cubelist = cubes.concatenate()
                print cubelist
                if len(cubelist[0].cell_methods) == 0:
                    cubelist[0].cell_methods = cubelist[1].cell_methods
                cubes = cubelist
            if 'GPCC' in infile[0]:
                method_con = iris.Constraint(cube_func = lambda c: \
                                any(m.method=='mean' for m in c.cell_methods))
                cubes = cubes.extract('precipitation_flux' & method_con)
            self.original_data['ppn'] = cubes.concatenate_cube()
            #if 'parent_source_id' in self.original_data['ppn'].attributes:
            #    if 'EC-Earth3' in self.original_data['ppn'].attributes['parent_source_id']:
            #        redefine_spatial_coords(self.original_data['ppn'])
            if 'model_id' in self.original_data['ppn'].attributes:
                 if ('ALADIN' in self.original_data['ppn'].attributes['model_id']) |\
                    ('ALARO' in self.original_data['ppn'].attributes['model_id']):
                     handle_ALADIN_Lambert_grid(self.original_data['ppn'])
            if self.original_data['ppn'].coords('time'):
                self.start_year = self.original_data['ppn'].coord('time').cell(0).point.year
                self.end_year = self.original_data['ppn'].coord('time').cell(-1).point.year
        # if ('alps' in infile[0].lower()) & (('n512' not in infile[0].lower()) | ('mi-ao438' not in infile[0].lower()) | (('12km' not in infile[0].lower()))):
        #    add_laea.load_laea(self.original_data['ppn'], '/project/hires_rcm/hadek/ALPS/RapdD_al05.etrs.laea_1989.nc')

        if maskfile: 
            self.original_data['mask'] = iris.load_cube(maskfile)
            if not isinstance(self.original_data['mask'].data, np.ma.MaskedArray):
                masked_cube = self.original_data['mask'].copy()
                masked_cube.data = np.ma.MaskedArray(self.original_data['mask'].data, mask=False)
                self.original_data['mask'] = masked_cube
            # print maskfile
            print 'mask.shape: ', self.original_data['mask'].shape
            print 'data.shape: ', self.original_data['ppn'].shape
        else: 
            if len(self.config['lsmasks'][res]) == 0:
                print 'generating false mask'
                filenames = glob.glob(self.config['infiles'])
                sample_cube = iris.load(filenames[0])
                latitude = sample_cube[0].coord('latitude')
                longitude = sample_cube[0].coord('longitude')
                cube_mask = Cube(np.ones(sample_cube[0].shape[1:]), dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
                self.original_data['mask'] = cube_mask
            else:
                print 'loading mask from config[lsmasks] defined in process_data.py'
                self.original_data['mask'] = iris.load_cube(self.config['lsmasks'][res])
        if self.original_data['ppn'].units != 1:
            self.original_data['ppn'] = check_precip_mm_or_kg_m2_s1(self.original_data['ppn'])


    def load_data(self, regions=None, season=None):
        '''
        populates self.hists ad self.totnum from JSON files for the regions in self or provided
        in the function input, for resolution in self and for files in self.config['histfile'] (different versions of same resolution)
        This loaded data will be used for plotting.

        :param list regions: regions to be loaded
        :param str season: season to be loaded
        '''

        metadata = {'totnum': self.totnum,
                    'nb_spatial_points': self.nb_spatial_points,
                    'start_year': self.start_year,
                    'end_year': self.end_year,
                    'centiles': self.centiles,
                    'yearly_time_series': self.yearly_time_series}

        for res in self.resolutions:
            self.hists[res] = {}
            self.hist2d[res] = {}
            self.totnum[res] = {}
            self.nb_spatial_points[res] = {}
            self.start_year[res] = {}
            self.end_year[res] = {}
            self.centiles[res] = {}
            self.yearly_time_series[res] = {}
            if not regions:
                regions = self.regions.keys()
            for r in regions:
                self.hists[res][r] = {}
                self.hist2d[res][r] = {}
                self.totnum[res][r] = {}
                self.nb_spatial_points[res][r] = {}
                self.start_year[res][r] = {}
                self.end_year[res][r] = {}
                self.centiles[res][r] = {}
                self.yearly_time_series[res][r] = {}
                for hfile in self.config['histfile'][res]:
                    print("{}_{}_{}_{}".format( res, r, season, hfile))
                    with open("{}{}_{}_{}_{}".format(self.config['jsondir'], res, r, season, hfile), 'r') as fh:
                        tmp = json.load(fh)
                        numpy_walk(tmp)
                        self.hists[res][r][hfile] = tmp['hists']
                        if 'hist3d' in tmp:
                            self.hist2d[res][r][hfile] = tmp['hist3d']
                        self.bins = tmp['bins']
                        if tmp['duration_bins'] is not None:
                            self.duration_bins = tmp['duration_bins']
                        for json_key in metadata.keys():
                            if json_key in tmp.keys():
                                if tmp[json_key] is not None:
                                    metadata[json_key][res][r][hfile] = tmp[json_key]

    def dump_data(self, targfile):
        '''
        Puts the self.hists, self.bins, self.duration_bins and self.totnum into a JSON file
        
        :param str targfile: File name of the target file that the histogram and bin data will be moved into 
        '''

        with open(targfile, 'w') as fh:
            output = {'bins': self.bins, 'hists': self.hists,
                      'hist3d': self.hist2d,
                      'duration_bins': self.duration_bins,
                      'totnum': self.totnum,
                      'nb_spatial_points': self.nb_spatial_points,
                      'start_year': self.start_year,
                      'end_year': self.end_year,
                      'centiles': self.centiles,
                      'yearly_time_series': self.yearly_time_series}
            json.dump(output, fh, indent=2, sort_keys=True, cls=NumpyAwareJSONEncoder)

    def model_averages(self):
        ''' 
        Calculates the average of the ensemble members to help with later statistical analysis
        
        '''
        model_av = {}
        for res in self.resolutions:
            model_av[res] = {}
            for reg in self.regions.keys():
                mav = []
                for i in range(len(self.hists[res][reg][self.config['histfile'][res][0]])):
                    step = []
                    for hf in self.config['histfile'][res]:
                        step.append(self.hists[res][reg][hf][i])
                    mav.append(np.mean(step))
                self.hists[res][reg]['average'] = np.array(mav)
            self.config['histfile'][res].append('average')

    def geostatial_limits(self):
        '''
        Calculates the geostatial longitude and latitude ranges for cubes with and without them in their metadata
        
        :param dict original_data: Contains the original mask and ppn that were loaded
        :return list geo_limits: List of geostatial limits extracted from cube
        '''
        if 'geostatial_lon_max' in self.original_data['ppn'].attributes.keys():
            geo_lon_max = self.original_data['ppn'].attributes['geostatial_lon_max']
        else:
            geo_lon_max = np.round(self.original_data['ppn'].coord('longitude').points.max())
        if 'geostatial_lon_min' in self.original_data['ppn'].attributes.keys():
            geo_lon_min = self.original_data['ppn'].attributes['geostatial_lon_min']
        else:
            geo_lon_min = np.round(self.original_data['ppn'].coord('longitude').points.min())
        if 'geostatial_lat_max' in self.original_data['ppn'].attributes.keys():
            geo_lat_max = self.original_data['ppn'].attributes['geostatial_lat_max']
        else:
            geo_lat_max = np.round(self.original_data['ppn'].coord('latitude').points.max())
        if 'geostatial_lat_min' in self.original_data['ppn'].attributes.keys():
            geo_lat_min = self.original_data['ppn'].attributes['geostatial_lat_min']
        else:
            geo_lat_min = np.round(self.original_data['ppn'].coord('latitude').points.min())
        geo_limits = [geo_lat_max, geo_lat_min, geo_lon_max, geo_lon_min]
        return geo_limits

    def extract_mask_ppn(self, reg, extension, ppn=None, mask=None, mask_or_not = True):
        '''
        Extracts the correct area from the mask and ppn given the region wanted
        
        :param str reg: Region to be selected
        :param int extension: How much extra to select beyond region for regridding purposes
        :param iris.cube ppn: Data to have region extracted from, default None means original data provided is used
        :param iris.cube mask: Mask to have region extracted from, default None means original mask provided is used
        :return iris.cube r_ppn: Regional data set
        :return iris.cube r_mask: Regional land-sea mask
        '''


        b = list(self.regions[reg].bounds)

        geo_limits = self.geostatial_limits()

        if self.original_data['ppn'].coord(
                'longitude').points.min() >= 0:  # Warning: this may cause bugs in regional configurations where the minimum latitude is above 0 but in a -180:180 system.
            print "assuming longitude in a 0->360 degree system"
            diff = 180
        else:
            diff = 0
        reg_lon_ce = iris.coords.CoordExtent('longitude', b[0] - extension, b[2] + extension)
        reg_lat_ce = iris.coords.CoordExtent('latitude', b[1] - extension, b[3] + extension)

        ### ensure that the region wanted is within the ppn bounds ###
        if reg_lat_ce.maximum <= geo_limits[0] and reg_lat_ce.minimum >= geo_limits[1] - 1 and reg_lon_ce.maximum <= \
                geo_limits[2] - diff + 3. and reg_lon_ce.minimum >= geo_limits[3] - diff - 1:
            if ppn is None:
                r_ppn = extract_region(self.original_data['ppn'], reg_lon_ce, reg_lat_ce)
            else:
                r_ppn = extract_region(ppn, reg_lon_ce, reg_lat_ce)
            if mask is None:
                r_mask = extract_region(self.original_data['mask'], reg_lon_ce, reg_lat_ce)
            else:
                r_mask = extract_region(mask, reg_lon_ce, reg_lat_ce)
            if not mask_or_not:
                print('Warning: no land-sea mask used, all points used')
                r_mask1 = r_mask.copy()
                if isinstance(r_mask.data, np.ma.MaskedArray):
                    r_mask1.data.mask = False
                else:
                    r_mask1.data[r_mask.data < 1] = 1
                r_mask = r_mask1
        else:
            print "region lon: {} to {}, lat: {} to {} is not in the cube".format(b[0], b[2], b[1], b[3])
            r_ppn = None
            r_mask = None

        return r_ppn, r_mask

    def gendata(self, regrid_res, reg, wet_or_dry_spells=None, thrs=None, hr2day=False,
                dry_window_allowed_within_wet_spell=None, testing=False, mask_or_not=True):
        '''
        Generates histogram and bin data for use in GUI
        
        :param str reg: Defines the region for which data will be generated
        :param str regrid_res: The resolution to regrid the data to
        :param str season: The season you wish to look at
        '''

        def lon(l, lon_range=None):
            ''' Adjusts longitudinal co-oridnate if falls outside normal range '''
            if lon_range is None: lon_range = (0., 360.)
            if l < lon_range[0]: l += 360
            if l > lon_range[1]: l -= 360
            return l

        def bound_cube(ppn_data):
            coord_names = [coord.name() for coord in ppn_data.dim_coords]
            try:
               lonname, = filter(lambda x: re.search(r'[Ll]on|x', x), coord_names)
            except:
               add_CORDEX_grid(ppn_data)
               coord_names = [coord.name() for coord in ppn_data.dim_coords]
               lonname, = filter(lambda x: re.search(r'[Ll]on|x', x), coord_names)
            latname, = filter(lambda x: re.search(r'[Ll]at|y', x), coord_names)
            '''Creates bounds for an unbounded '''
            if not ppn_data.coord(latname).has_bounds():
                lat_coord = ppn_data.coord(latname)
                lat_coord.guess_bounds()
            if not ppn_data.coord(lonname).has_bounds():
                lon_coord = ppn_data.coord(lonname)
                lon_coord.guess_bounds()

        def regrid(regrid_res, reg, r_ppn, ldmask, mask_or_not):
            ''' Regrids the mask and ppn on a given resolution '''
            print 'Regridding at this resolution: ', regrid_res

            if r_ppn is not None:
                config1 = self.config.copy()
                config1['infiles'] = self.config['lsmasks'][regrid_res]
                Newgrid_class = PrecipExtremes(srex_reg, config1, regrid_res)
                Newgrid_class.load_mask_ppn(regrid_res)
                new_ppn, r_mask = Newgrid_class.extract_mask_ppn(reg, 3, mask_or_not=mask_or_not)
                if isinstance(r_ppn, iris.cube.Cube):
                    print ldmask.shape, r_ppn.shape
                #      landmask = iris.util.broadcast_to_shape(ldmask.data, r_ppn.shape, (1, 2))
                #      r_ppn = mask_cube_where(r_ppn, landmask < 1)
                if new_ppn is None:
                    print 'Region ' + reg + ' not in regridding area. Skipping'
                    r_ppn = None
                else:
                    bound_cube(new_ppn)
                    if r_ppn.coord_system() == None:
                        r_ppn.coord('latitude').coord_system = new_ppn.coord_system()
                        r_ppn.coord('longitude').coord_system = new_ppn.coord_system()
                    r_ppn02 = r_ppn.regrid(new_ppn, iana.AreaWeighted(mdtol=0.2))
                    r_ppn = r_ppn02
                    if not isinstance(r_ppn.data, np.ma.core.MaskedArray):
                        r_ppn.data = np.ma.MaskedArray(r_ppn.data, mask=False)
                r_ppn, r_mask = self.extract_mask_ppn(reg, 0, ppn=r_ppn, mask=r_mask, mask_or_not=mask_or_not)

                if np.sum(r_mask.shape) != np.sum(r_ppn[0, :, :].shape):
                    if r_mask.coord_system() == None:
                        r_mask.coord('latitude').coord_system = r_ppn.coord_system()
                        r_mask.coord('longitude').coord_system = r_ppn.coord_system()
                    r_mask.regrid(r_ppn[0, :, :], iana.AreaWeighted())  # r_mask[:,:r_ppn[0,0,:].shape[0]]
                    r_mask = r_mask.intersection(
                        longitude=(r_ppn.coord('longitude').cell(0).point, r_ppn.coord('longitude').cell(-1).point),
                        latitude=(r_ppn.coord('latitude').cell(0).point, r_ppn.coord('latitude').cell(-1).point))

            else:
                r_mask = None
            #            i = 100
            #            plt.figure()
            #            ax = plt.subplot(1,3,1, projection = ccrs.PlateCarree())
            #            iplt.contourf(orig_r_ppn[i,:,:])
            #            ax.coastlines()
            #            plt.title('Non-regridded')
            #            ax1 = plt.subplot(1, 3, 2, projection = ccrs.PlateCarree())
            #            iplt.contourf(r_ppn_old[i,:,:])
            #            ax1.coastlines()
            #            plt.title('Original method')
            #            ax1 = plt.subplot(1, 3, 3, projection = ccrs.PlateCarree())
            #            iplt.contourf(r_ppn[i,:,:])
            #            ax1.coastlines()
            #            plt.title('New method')
            #            plt.show()
            #            plt.figure()
            #            ax = plt.subplot(1,1,1, projection = ccrs.PlateCarree())
            #            iplt.contourf(orig_r_ppn[100,:,:], levels = np.arange(0,200, 25))
            #            ax.coastlines()
            #            plt.colorbar()
            #            plt.title(self.resolutions[0])
            #            plt.show()

            return r_ppn, r_mask

        print 'original cube:', self.original_data['ppn']
        bound_cube(self.original_data['ppn'])
        r_ppn = self.original_data['ppn']
        r_mask = self.original_data['mask']

        if reg != 'masked':  # reg = 'masked means that the mask for the region will come from r_mask directly'
          if 'latitude' in r_mask.coords():
            if regrid_res:
                print 'regridding to:', regrid_res
                r_ppn, r_mask = self.extract_mask_ppn(reg, 3, mask_or_not=mask_or_not)
            else:
                r_ppn, r_mask = self.extract_mask_ppn(reg, 0., mask_or_not=mask_or_not)
        print 'extracted cube to be used for regional histograms:', r_ppn
        print 'cube starting ', r_ppn.coord('time').cell(0)
        print 'cube ending ', r_ppn.coord('time').cell(-1)

        if r_ppn is None:
            print 'Region ' + reg + ' not in data provided. Skipping'
            self.hists = np.zeros(self.bins.shape)
            # return self.hists
        else:
            if hr2day:  # tranforms the hourly cube to a daily cube
                print 'aggregating hourly to daily...'
                daily_cubes = iris.cube.CubeList()
                for year in range(r_ppn.coord('time').cell(0).point.year, r_ppn.coord('time').cell(-1).point.year + 1):
                    time_con1 = iris.Constraint(time=lambda cell: cell.point.year == year)
                    subcub = r_ppn.extract(time_con1)
                    daily_cubes.append(hourly2daily(subcub))
                r_ppn = daily_cubes.concatenate_cube()
            if regrid_res:
                r_ppn, r_mask = regrid(regrid_res, reg, r_ppn, r_mask, mask_or_not)
            if r_ppn is None:
                print 'Region or season not in data provided. Skipping'
                self.hists = np.zeros(self.bins.shape)
                # return self.hists

            # sometimes the r_mask shape is slightly different from the r_ppn shape
            rr = int(random.random() * 100)
            if np.sum(r_mask.shape) != np.sum(r_ppn[0, :, :].shape):
                print 'regridding mask from mask shape {} to r_ppn shape {}'.format(r_mask.shape, r_ppn.shape)
                lat_coord, lon_coord = cube_primary_xy_coord_names(r_mask)
                if r_mask.coord_system() == None:
                    r_mask.coord(lon_coord).coord_system = r_ppn.coord_system()
                    r_mask.coord(lat_coord).coord_system = r_ppn.coord_system()
                if not r_mask.coord(lat_coord).has_bounds():
                    r_mask.coord(lon_coord).guess_bounds()
                    r_mask.coord(lat_coord).guess_bounds()
                lat_coord, lon_coord = cube_primary_xy_coord_names(r_ppn)
                if not r_ppn.coord(lon_coord).has_bounds():
                    r_ppn.coord(lon_coord).guess_bounds()
                    r_ppn.coord(lat_coord).guess_bounds()
                if r_ppn.coord_system() == None:
                    if lon_coord == 'longitude':
                       r_ppn.coord(lon_coord).coord_system = iris.coord_systems.GeogCS(6371229.0)
                       r_ppn.coord(lat_coord).coord_system = iris.coord_systems.GeogCS(6371229.0)
                       print 'assuming r_ppn coord system is iris.coord_systems.GeogCS(6371229.0)'
                # qplt.contourf(r_ppn[0, :, :])
               
                if reg == 'masked': 
                 if 'model_id' in r_ppn.attributes:
                   if ('ALADIN' in r_ppn.attributes['model_id']) | ('ALARO' in r_ppn.attributes['model_id']):
                       r_ppn.coord('projection_x_coordinate').convert_units('m')
                       r_ppn.coord('projection_y_coordinate').convert_units('m')
                   if ('NUIM-WRF341E' in r_ppn.attributes['model_id']) | ('IPSL-INERIS-WRF331F' in r_ppn.attributes['model_id']):
                       r_ppn.remove_coord('grid_latitude')        
                       r_ppn.remove_coord('grid_longitude')
                       add_CORDEX_grid(r_ppn)
                   if 'ICTP-RegCM4-3' in r_ppn.attributes['model_id']:
                       el_axis = r_ppn.coord('projection_x_coordinate').coord_system.ellipsoid.semi_major_axis
                       r_ppn.coord('projection_x_coordinate').coord_system.ellipsoid.semi_minor_axis = el_axis
                       r_ppn.coord('projection_y_coordinate').coord_system.ellipsoid.semi_minor_axis = el_axis
                r_mask = r_mask.regrid(r_ppn[0, :, :], iana.Linear('mask'))  # r_mask[:,:r_ppn[0,0,:].shape[0]]
                   
                #if reg is not 'masked':
                #    fig = plt.figure()
                #    axx = fig.add_subplot(111, projection=ccrs.PlateCarree())
                #    qplt.contourf(r_ppn[0, :, :])
                #    axx.coastlines()
                #    plt.savefig(
                #        '/home/users/sberthou/PrecipExtremes/images/rppn_{}_{}.png'.format(reg,
                #                                                                                                rr))

            regenmask = False
            if not isinstance(r_mask.data, np.ma.core.MaskedArray):
                regenmask = True
            elif r_mask.data.mask.shape == ():
                regenmask = True
            if regenmask:
                print 'r_mask not masked array'
                if r_mask.data.max() == 1.0:
                    print 'generating ocean mask from land-sea fraction'
                    r_mask = mask_cube_where(r_mask, r_mask.data < 0.1)
                    if len(r_mask.data.mask.shape) == 0:
                        r_mask.data = np.ma.MaskedArray(r_mask.data, mask=False)
                else:
                    r_mask1 = r_mask.copy()
                    r_mask1.data = np.ma.MaskedArray(r_mask.data, mask=False)
                    r_mask = r_mask1

            if reg != 'masked': # if using Polygon defined regions, creating the mask for it
                print 'r_mask', r_mask
                if r_mask.coords('grid_latitude'):
                    rot_pole = r_mask.coord('grid_latitude').coord_system.as_cartopy_crs()
                    ll = ccrs.Geodetic()
                for i in range(r_mask.data.shape[0]):
                    for j in range(r_mask.data.shape[1]):
                        if r_mask.coords('grid_latitude'):
                           rlon = r_mask.coord('grid_longitude').points[j]
                           rlat = r_mask.coord('grid_latitude').points[i] 
                           lon_lat = ll.transform_point(rlon, rlat, rot_pole)
                           m_point = Point(lon_lat[0], lon_lat[1])
                        else:
                           m_point = Point(lon(r_mask.coord('longitude').points[j], lon_range=(-180., 180.)),
                                        r_mask.coord('latitude').points[i],)
                        if not(self.regions[reg].contains(m_point)):
                            if not r_mask.data.mask[i, j]:
                                r_mask.data.mask[i, j] = True

            bound_cube(r_mask)
#            ### uncommnent the following lines if you want to check that the mask is well defined
#            fig = plt.figure()
#            axx = fig.add_subplot(111, projection=ccrs.PlateCarree())
#            qplt.contourf(r_mask)
#            axx.coastlines("50m")
#            plt.savefig('/home/users/sberthou/PrecipExtremes/images/{}_{}.png'.format(reg, rr))
            print 'number of non masked points:', np.ma.count(r_mask.data)
            self.nb_spatial_points = np.ma.count(r_mask.data)
            area_weights = np.ma.MaskedArray(np.ones(
                    r_mask.shape))  # area weights doesn't work for non regular grids; these grids usually are area equal, or almost area equal.
            area_weights.mask = r_mask.data.mask

            self.hists = {}

            ########### change the bins (in kg/m2/s1) into the file units ##############
            if r_ppn.units in ['mm/day', 'mm day-1']:
                bins1 = self.bins.copy()
                self.bins = self.bins * 86400.
            elif r_ppn.units in ['mm/h', 'mm h-1', 'mm (h)-1', 'mm/hour']:
                bins1 = self.bins.copy()
                self.bins = self.bins * 3600.
            elif r_ppn.units in ['mm (3h)-1', 'mm/3h']:
                bins1 = self.bins.copy()
                self.bins = self.bins * 3600. * 3.
            elif r_ppn.units in ['mm']:
                print 'warning, if you are not using AMMA-CATCH 15mn, this will be wrong'
                bins1 = self.bins.copy()
                self.bins = self.bins * 15. * 60.
            elif r_ppn.units in ['kg m-2 s-1']:
                pass
            else:
                raise TypeError('r_ppn units not in mm/day, mm/h, mm/3h or mm')
            print r_ppn.units, self.bins

            if wet_or_dry_spells is None:
                self.gen_hist1d(r_ppn, area_weights)
                lat_name, lon_name = cube_primary_xy_coord_names(r_ppn)
                self.centiles_over_masked_region(r_ppn, area_weights, lat_name, lon_name)
                self.annual_time_series(r_ppn, area_weights, lat_name, lon_name)
            else:
                self.amount_duration_hist(r_ppn, area_weights, wet_or_dry_spells, thrs,
                                          dry_window_allowed_within_wet_spell, testing=testing)
            ########### put the bins back in kg/m2/s1 for saving ##############
            if r_ppn.units in ['mm/day', 'mm', 'mm day-1', 'mm/h', 'mm h-1', 'mm (h)-1', 'mm/hour',
                               'mm (3h)-1']:  # , 'kg m-2 s-1']: #just for Giorgia's runs
                self.bins = bins1
            # return self.hists, area_weights

    def gen_hist1d(self, r_ppn, area_weights):
        yr = 0
        step_count = 1
        yrange = range(r_ppn.coord('time').cell(0).point.year, r_ppn.coord('time').cell(-1).point.year+1)
        self.hists = {}
        for year in yrange:
            self.hists[year] = [0] * (len(self.bins) - 1)
            print year
            time_con = iris.Constraint(time=lambda cell: cell.point.year == year)
            cube = r_ppn.extract(time_con)
            for s in cube.slices_over(['time']):  # slices(['latitude','longitude']):
                s_data_masked = mask_cube_where(s, area_weights.mask)
                area_weight_resize = area_weights.compressed()
                s_data_resize = s_data_masked.data.compressed()
                if len(s_data_resize) < len(area_weight_resize):
                    # print 'suspicious: len(s.data.compressed)', len(s_data_resize), \
                    #      'not equal to len(area_weights.compressed())', len(area_weight_resize)
                    area_weights1 = area_weights.copy()
                    area_weights1.mask = s_data_masked.data.mask
                    area_weight_resize = area_weights1.compressed()
                elif len(s_data_resize) > len(area_weight_resize):
                    # print 'suspicious: len(s.data.compressed)', len(s_data_resize), \
                    #		  'not equal to len(area_weights.compressed())', len(area_weight_resize)
                    s_data_resize = s_data_resize[0:len(area_weight_resize)]

                c, _ = np.histogram(s_data_resize, self.bins, weights=area_weight_resize)
                self.hists[year] += c
                step_count += 1

    def centiles_over_masked_region(self, r_ppn, area_weights, lat_name='latitude', lon_name = 'longitude'):
        centiles = {}

        r_ppn_masked = mask_3Dcube_with_2Dmask(r_ppn, area_weights.mask)
        for percentile in [60, 65, 66, 67, 68, 70, 80, 85, 88, 90, 95, 99, 99.5, 99.9, 99.99, 99.999]:
            centiles['p' + str(percentile)] = np.float(r_ppn_masked.collapsed(('time', lat_name, lon_name), \
                                                                              iana.PERCENTILE,
                                                                              percent=percentile).data)  # , weights=area_weights).data)
        centiles['units'] = str(r_ppn.units)
        self.centiles = centiles

    def annual_time_series(self, r_ppn, area_weights, lat_name='latitude', lon_name='longitude'):
        r_ppn_masked = mask_3Dcube_with_2Dmask(r_ppn, area_weights.mask)
        icat.add_year(r_ppn_masked, 'time')
        yearly_sum = r_ppn_masked.aggregated_by('year', iana.SUM)
        yearly_time_series = yearly_sum.collapsed((lat_name, lon_name), iana.MEAN)
        self.yearly_time_series = yearly_time_series.data

    def amount_duration_hist(self, cube, area_weights, wet_or_dry_spells, thrs,
                             dry_window_allowed_within_wet_spell=None, testing=False):
        '''
        Uses the function notmasked_contiguous on a masked array (below or above config['thrs']) 
        and returns bins of intensity and duration of wet or dry spells.

        :param iris.cube cube: 3D data (time, lat, lon)
        :param np.array area_weights: 
        :param str wet_or_dry_spells: 'dry' or 'wet'
        :param float thrs: in mm/day
        '''

        cube_masked = cube.copy()
        # if testing:
        #     fig = plt.figure()
        #     plt.contourf(area_weights)
        #     plt.title('area mask used')
        #     plt.savefig('/data/users/sberthou/PrecipExtremes/srex_prext_other_regions/images/area_weight_{:d}.png'.format(int(np.random.rand()*100)))
        # print area_weights, area_weights.shape
        area_weight1 = iris.util.broadcast_to_shape(area_weights, cube.shape, (1, 2))
        cube_masked.data = np.ma.MaskedArray(cube.data, mask=area_weight1.mask)
        # if testing:
        #     fig = plt.figure()
        #     iplt.contourf(cube1[0, :, :])
        #     plt.title('area mask used')
        #     plt.savefig('/data/users/sberthou/PrecipExtremes/srex_prext_other_regions/images/cube_{:d}.png'.format(int(np.random.rand()*100)))
        if testing:
            print('nb non masked data for region:', len(cube_masked.data.compressed()))
        if cube_masked.units in ['mm/h', 'mm h-1', 'mm (h)-1', 'mm/hour']:
            thrs /= 24.
            print 'thrs = {}mm/h'.format(thrs)
        elif cube_masked.units in ['mm (3h)-1']:
            thrs /= 8.
            print 'thrs = {}mm/3h'.format(thrs)
        elif cube_masked.units in ['kg/m2/s', 'kg m-2 s-1']:
            thrs /= 86400.
            print 'thrs = {}kg/m2/s'.format(thrs)
        elif cube_masked.units in ['mm']:
            thrs /= 96.
            print 'thrs = {}mm/15mn'.format(thrs)
        else:
            print 'thrs = {}mm/day'.format(thrs)
        if wet_or_dry_spells == 'dry':
            cube_masked.data = np.ma.masked_where(cube.data > thrs, cube_masked.data)
        else:
            cube_masked.data = np.ma.masked_where(cube.data < thrs, cube_masked.data)

        yrange = range(cube.coord('time').cell(0).point.year, cube.coord('time').cell(-1).point.year + 1)
        self.hist2d = {}
        self.hist = {}
        for year in yrange:
            print year
            time_con = iris.Constraint(time=lambda cell: cell.point.year == year)
            cube_year = cube.extract(time_con)
            cube_masked_year = cube_masked.extract(time_con)

            if testing:
                print('after masking with thrs: ', len(cube_masked.data.compressed()))
                cube2 = cube[:, 0:2, 0]
                if wet_or_dry_spells == 'wet':
                    series_to_test = np.array(
                        (0, thrs * 5, thrs * 0.1, thrs * 5, 0, 0, thrs * 5, thrs * 0.1, thrs * 5, 0, 0,
                         thrs * 0.5, thrs * 0.1, thrs * 4))
                else:
                    series_to_test = np.array(
                        (0, thrs * 5, thrs * 0.1, thrs * 2, 0, 0, thrs * 5, thrs * 0.1, thrs * 5, 0, 0,
                         thrs * 0.5, thrs * 0.1, 0, 0, 0, 0, thrs * 2))
                cube2.data[:, 0] = np.append(series_to_test, np.zeros(cube.shape[0] - len(series_to_test)))
                cube2.data[:, 1] = np.append(series_to_test, np.zeros(cube.shape[0] - len(series_to_test)))
                cube_masked = cube2.copy()
                self.bins = [0, thrs * 5, 2000]
                self.duration_bins = [0, 2, 4, cube.shape[0]]
                self.drought_bins = [0, 2, 4, 8, cube.shape[0]]
                if wet_or_dry_spells == 'dry':
                    cube_masked.data = np.ma.masked_where(cube2.data > thrs, cube_masked.data)
                else:
                    cube_masked.data = np.ma.masked_where(cube2.data < thrs, cube_masked.data)
                print cube_masked.data[0:10]
                tmp = cube_masked.data
                tmp_not_masked = cube2.data
            else:
                tmp_not_masked = cube_year.data.reshape(cube_year.data.shape[0],
                                                        cube_year.data.shape[1] * cube_year.data.shape[2])
                tmp = cube_masked_year.data.reshape(cube_masked_year.data.shape[0],
                                                    cube_masked_year.data.shape[1] * cube_masked_year.data.shape[2])
                # if you want to add area weights:
                weights1d = area_weights.reshape(
                    cube_masked_year.data.shape[1] * cube_masked_year.data.shape[2]).filled(
                    0)  # give 0 shape where landmask
            slices = np.ma.notmasked_contiguous(tmp, axis=0)
            if testing:
                print(slices)
            # hist2d[i] = np.zeros(cube1.data.shape[1]*cube1.data.shape[2])
            durations = np.ones(0)
            intensities = np.ones(0)
            following_drought_length = np.ones(0)
            count = 0
            for i, sli in enumerate(slices):
                # print sli
                if isinstance(sli, list):
                    if len(sli) > 1:
                        ##################### change mean to max for peak intensity #######################
                        # print np.array([s.stop-s.start for s in sli])
                        # print np.array([tmp[s].mean() for s in sli])
                        if dry_window_allowed_within_wet_spell:
                            critic_intens = thrs * 3
                            num = 0
                            while num < len(sli):
                                num_start = num
                                fol_drought_len = dry_window_allowed_within_wet_spell
                                fol_wet_intensity = critic_intens
                                cdt_next_wet_spell_large_enough = True
                                while (
                                        fol_drought_len <= dry_window_allowed_within_wet_spell) and cdt_next_wet_spell_large_enough and (
                                        num < len(sli)):
                                    dura = sli[num].stop - sli[num_start].start
                                    inten = tmp[sli[num_start].start:sli[num].stop, i].sum()
                                    if num == len(sli) - 1:
                                        fol_drought_len = tmp.shape[0] - sli[num].stop
                                        fol_wet_intensity = 0
                                    else:
                                        fol_drought_len = sli[num + 1].start - sli[num].stop
                                        fol_wet_intensity = tmp_not_masked[sli[num + 1].start - 1, i]
                                        # print 'num+1', num+1, sli[num+1].start-1, tmp_not_masked[sli[num+1].start-1, i]
                                    if wet_or_dry_spells == 'dry':
                                        cdt_next_wet_spell_large_enough = (fol_wet_intensity <= critic_intens)
                                        # print cdt_next_wet_spell_large_enough, fol_wet_intensity, critic_intens
                                    num = num + 1
                                durations = np.append(durations, dura)
                                intensities = np.append(intensities, inten)
                                following_drought_length = np.append(following_drought_length, fol_drought_len)
                        else:
                            durations = np.append(durations, np.array([s.stop - s.start for s in sli]))
                            intensities = np.append(intensities, np.array([tmp[s.start:s.stop, i].sum() for s in sli]))
                            following_drought_length = np.append(following_drought_length,
                                                                 np.array([sp1.start - s.stop for s, sp1 in
                                                                           zip(sli[:-1], sli[1:])]))
                            following_drought_length = np.append(following_drought_length,
                                                                 cube_masked_year.data.shape[0] - sli[-1].stop)
                    else:
                        #print i, sli
                        if len(sli) == 1:
                            if tmp[sli[0].start:sli[0].stop, i].sum() is not np.ma.masked:
                                durations = np.append(durations, sli[0].stop - sli[0].start)
                                intensities = np.append(intensities, tmp[sli[0].start:sli[0].stop, i].sum())
                                following_drought_length = np.append(following_drought_length,
                                                                     cube_masked_year.data.shape[0] - sli[0].stop)
                                # may add weights according to area of each cell here:
                                # np.append(weights, weights1d[i]*np.ones(len(sli)) # something like that, put 0 weight where ocean.
                        else:
                            if sli == []:
                                count += 1
                            else:
                                if tmp[sli.start:sli.stop, i].sum() is not np.ma.masked:
                                    durations = np.append(durations, sli.stop - sli.start)
                                    intensities = np.append(intensities, tmp[sli.start:sli.stop, i].sum())
                                    following_drought_length = np.append(following_drought_length,
                                                                         cube_masked_year.data.shape[0] - sli.stop)

            self.totnum[year] = len(durations)
            print 'totnum', self.totnum
            print 'number of empty slices', count
            print len(durations)
            hist2d, dur_bins, intens_bins = np.histogram2d(durations, intensities,
                                                           bins=(self.duration_bins, self.bins),
                                                           normed=False)
            self.hists[year] = hist2d
            if wet_or_dry_spells in ['wet', 'dry']:
                hist3d, edges = np.histogramdd((durations, intensities, following_drought_length),
                                               bins=(self.duration_bins, self.bins, self.drought_bins),
                                               normed=False)
                self.hist2d[year] = hist3d
                self.edges = edges
                if testing:
                    print(hist3d.shape)
                    for i in range(hist3d.shape[0]):
                        print(i, 'duration')
                        print(hist3d[i, :, :])
                    print(edges[0])
                    print(edges[1])
                    print(edges[2])
                    raise UserWarning('stop because testing')


    def save_subplot(self, reg_todo, plot_type, x_lim, colorlwdict, figdir, average=False, season=None, bigregion=None,
                     subregion_files=None, nb_freq_per_day=1., diff_btw_res=None, dict_names={}):
        '''
        Save all the plots

        :param dict reg_todo: Contains all the region shapes
        :param str plot_type: Type of data being plotted to be included in filename: contrib, frac, normed, normed_wet, centiles
        :param int x_lim: Sets the max of the x-axis
        :param dict colorlwdict: contains [ line color, linewidth, linetype] for each model
        :param str figdir: directory where the figure will be saved
        :param bool average: not used anymore, used to be to take the average of different ensemble members
        :param str season: e.g 'jja'
        :param dict subregion_files: for interface_plot_european_masked_subregions: uses masks to define regions
        :param float nb_freq_per_day: 24. if hourly, 8. if 3 hourly, 1. if daily
        :param list diff_btw_res: diff_btw_res = ['dataset_ref', 'dataset1', 'dataset_ref', 'dataset2',]
                                  will plot dataset1 - dataset_ref and dataset2 - dataset_ref
        :param dict dict_names: optional, if you want a different label name than the resolution name
        '''
        from matplotlib import style
        style.use('seaborn-colorblind')
        freqs = {'1': 'd', '8': '3h', '24': 'h', '96': '15mn'}
        print 'diff_btw_res', diff_btw_res, bool(diff_btw_res)
        if diff_btw_res:
            # plot the difference: diff_btw_res[1]-diff_btw_res[0], diff_btw_res[3]-diff_btw_res[2], ...
            res_todo = diff_btw_res[1::2]  # take the first and then every 2 datasets
            diff_res_todo = diff_btw_res[::2]  # take reference datasets every 2 datasets.
        else:
            res_todo = self.resolutions
        hfiter = False
        season_print = season
        if isinstance(season, list):
            season_print = ['_'.join(s for s in season)]
        if diff_btw_res:
            target = "{}/{}_{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season_print,
                                                int(24. / nb_freq_per_day), diff_btw_res[0])
        else:
            target = "{}/{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season_print,
                                             int(24. / nb_freq_per_day))
        regions = reg_todo
        if reg_todo == 'masked':
            hfiter = True
            reg_todo = self.hists[self.hists.keys()[0]]['masked'].keys()
            regions = None
            target = "{}/{}_{}_{}_{}_{}_{}".format(figdir, '_'.join(s[4:-5] for s in reg_todo), res_todo[0], plot_type,
                                                   nb_freq_per_day, freqs[str(int(nb_freq_per_day))], season_print)

        print target
        num_plots = len(reg_todo)
        # ### plots the regions of interest first ### does not plot the regions if only 1 region.
        fig, subpl, numcol, start_plot, num_lines, col_order = set_up_fig(num_plots, reg_todo=regions,
                                                                          bigregion=bigregion,
                                                                          subregion_files=subregion_files,
                                                                          season=season)

        if regions:
            if isinstance(bigregion, tuple):
                if subpl:
                    subpl.set_extent(bigregion, ccrs.PlateCarree())  # Africa
        for regcount, reg in enumerate(reg_todo):
            numplot = start_plot + regcount
            # new subplot for each region
            subpl2 = fig.add_subplot(num_lines, numcol, numplot)
            if season:
                season = reg
            print reg
            if reg is not None:
                bin_oc_tot = []
                arrow_dict = {}
                for nb_r, r in enumerate(res_todo):
                    if iris.util.approx_equal(self.bins[2] - self.bins[1], self.bins[21] - self.bins[20]): # traditional distrib with regular bins
                        bins1 = self.bins * 86400. / nb_freq_per_day
                    else:  # Martin distribution
                        print 'Martin distribution, assuming given in kg/m2/s, check in process_data.py', self.bins.shape
                        bins1 = self.bins * 86400. / nb_freq_per_day
                    if nb_r==0:
                        for xset in bins1:
                            plt.axvline(x=xset, linewidth=0.3, color='k')
                    if average:
                        #		runs=self.config['histfile'][r]['average']
                        hf = 'average'
                    else:
                        bin_mean = (bins1[:-1] + bins1[1:]) / 2.
                        ####### get the histograms to plot #########
                        if diff_btw_res:
                            hist_to_plot1, _, _, _ = self.get_hist2plot(r, reg, hfiter)
                            hist_to_plot1 = get_mean_10_year_change(hist_to_plot1[0])
                            print 'hist_to_plot1', hist_to_plot1.shape
                            hist_to_plotref, _, _, _ = self.get_hist2plot(diff_res_todo[nb_r], reg, hfiter)
                            hist_to_plotref = get_mean_10_year_change(hist_to_plotref[0])
                            print 'hist_to_plotref', hist_to_plotref.shape
                            if plot_type == 'normed':
                                hist_to_plot_fin = hist_to_plot1 / hist_to_plot1.sum() - hist_to_plotref / hist_to_plotref.sum()
                            elif plot_type == 'normed_wet':
                                if nb_freq_per_day == 1.:
                                    hist_to_plot_fin = (hist_to_plot1[1:] / hist_to_plot1[1:].sum() - hist_to_plotref[
                                                                                                      1:] / hist_to_plotref[
                                                                                                            1:].sum()) / (
                                                               hist_to_plotref[1:] / hist_to_plotref[
                                                                                     1:].sum()) * 100
                                    bins_to_plot = bin_mean[1:]
                                else:
                                    hist_to_plot_fin = (hist_to_plot1[2:] / hist_to_plot1[2:].sum() - hist_to_plotref[
                                                                                                      2:] / hist_to_plotref[
                                                                                                            2:].sum()) / (
                                                               hist_to_plotref[2:] / hist_to_plotref[
                                                                                     2:].sum()) * 100
                                    bins_to_plot = bin_mean[2:]

                            elif plot_type == 'frac':

                                hist_to_plot_fin = (hist_to_plot1 * bin_mean / sum(
                                    hist_to_plot1 * bin_mean) - hist_to_plotref * bin_mean / sum(
                                    hist_to_plotref * bin_mean)) * 100
                                hist_to_plot_fin_ref = hist_to_plotref * bin_mean / sum(
                                    hist_to_plotref * bin_mean) * 100
                            elif plot_type == 'contrib':
                                hist_to_plot_fin_ref = hist_to_plotref * bin_mean / hist_to_plotref.sum()
                                hist_to_plot_fin = hist_to_plot1 * bin_mean / hist_to_plot1.sum() - hist_to_plot_fin_ref
                                centiles = self.get_centiles(r, reg, hfiter)
                                centiles_ref = self.get_centiles(diff_res_todo[nb_r], reg, hfiter)
                            else:
                                hist_to_plot_fin = hist_to_plot1 - hist_to_plotref
                        else:
                            hist_to_plot, _, _, _ = self.get_hist2plot(r, reg, hfiter)
                            hist_to_plot = get_mean_10_year_change(hist_to_plot[0])
                            if plot_type == 'normed':
                                hist_to_plot_fin = hist_to_plot / hist_to_plot.sum()
                            elif plot_type == 'normed_wet':
                                if nb_freq_per_day == 1.:
                                    hist_to_plot_fin = hist_to_plot[1:] / hist_to_plot[1:].sum()
                                    bins_to_plot = bin_mean[1:]
                                else:
                                    hist_to_plot_fin = hist_to_plot[2:] / hist_to_plot[2:].sum()
                                    bins_to_plot = bin_mean[2:]
                            elif plot_type == 'frac':
                                hist_to_plot_fin = hist_to_plot * bin_mean / sum(hist_to_plot * bin_mean) * 100
                            elif plot_type == 'contrib':
                                hist_to_plot_fin = hist_to_plot * bin_mean / hist_to_plot.sum()
                            elif plot_type == 'centiles':
                                centiles = self.get_centiles(r, reg, hfiter)
                            else:
                                hist_to_plot_fin = hist_to_plot
                    ########### now plot, you may want to changes some things in here ##############
                    if plot_type == 'normed':
                        print 'total sum of frequencies for {}, {}: {}'.format(r, reg, hist_to_plot_fin.sum())
                        subpl2.plot(bin_mean[1:], hist_to_plot_fin[1:], '-',
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1],
                                    linestyle=colorlwdict[r][2],
                                    label='{}'.format(r))  # [:7]))
                        ylabel = 'Normed frequency'
                        if diff_btw_res:
                            #pass
                            subpl2.set_ylim([-0.02, 0.02])
                        else:
                            subpl2.set_ylim([0.000001, 0.2])  # daily Africa
                            subpl2.set_ylim([0.0000001, 0.1])  # daily Africa
                            # subpl2.set_ylim([0.00000002, 0.1])#[0.0000001,0.005])#[0.001, 0.5])#[0.0000001,0.005])#[0.0000001, 1])
                    elif plot_type == 'normed_wet':
                        print 'total sum of frequencies for {}, {}: {}'.format(r, reg, hist_to_plot_fin.sum())
                        subpl2.plot(bins_to_plot, hist_to_plot_fin, '-',
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1],
                                    linestyle=colorlwdict[r][2],
                                    label='{}'.format(r[:7]))
                        ylabel = 'Normed frequency (wet hours)'
                        if diff_btw_res:
                            subpl2.set_ylim([-500, 500])
                        else:
                            subpl2.set_ylim([0.000001, 1])
                    elif plot_type == 'frac':
                        subpl2.plot(bin_mean, hist_to_plot_fin,
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1], linestyle=colorlwdict[r][2],
                                    label='{}'.format(r))
                        subpl2.plot(bin_mean, hist_to_plot_fin_ref, color='k', label=r[:-7])
                        ylabel = 'Fractional contribution (\%) '
                        if diff_btw_res:
                            bin_oc_tot.append(sum(hist_to_plot1 * bin_mean) - sum(hist_to_plotref * bin_mean))
                            subpl2.set_ylim([-0.8, 0.8])
                        else:
                            subpl2.set_ylim([0, 3])
                            # plt.gca().set_yticks(range(4))
                            # plt.gca().set_yticklabels(range(4))
                    elif plot_type == 'contrib':
                        # if r is not 'TRMM':
                        lab = r
                        if r == 'N512r':
                            lab = 'R25'
                        if lab in dict_names.keys():
                            lab = dict_names[lab]
                        #if r is not 'TRMM':

                        if (diff_btw_res is not None) and (r == diff_res_todo[nb_r]): #in ['TRMM', 'CMORPH', 'CMORPH1']):
                            pass
                        else:
                            subpl2.plot(bin_mean, hist_to_plot_fin,
                                            color=colorlwdict[r][0], linewidth=colorlwdict[r][1],
                                            linestyle=colorlwdict[r][2], label='{}'.format(lab))
                            if r in ['SEA', 'AMZ']:
                                subpl2.set_ylim([0, 0.21])
                            else:
                                subpl2.set_ylim([0, 0.23]) #0.012]) #0.17]) #0.023])
                            if 24. / nb_freq_per_day == 1:
                                ylabel = 'Contribution (mm/h)'
                            elif nb_freq_per_day == 1:
                                ylabel = 'Contribution (mm/day)'
                            else:
                                ylabel = 'Contribution (mm/{}h)'.format(int(24. / nb_freq_per_day))

                        # subpl2.fill_between(bin_mean, hist_to_plot_fin, 0, color=colorlwdict[r][0], alpha=0.5)
                        # one_third, two_thirds = plot_one_third_two_thirds_vlines(subpl2, hist_to_plot_fin, bin_mean,
                        #                                      color=colorlwdict[r][0])
                        if (diff_btw_res is not None):
                            labref = diff_res_todo[nb_r]
                            if labref in dict_names.keys():
                                labref = dict_names[labref]
                            if (r != diff_res_todo[nb_r]):
                                subpl2.fill_between(bin_mean, hist_to_plot_fin, 0, color=colorlwdict[r][0], alpha=0.5)

                                one_third, two_thirds = plot_one_third_two_thirds_vlines(subpl2, hist_to_plot_fin_ref, bin_mean,
                                                             color=colorlwdict[labref][0])
                            # if r == 'N512r_future':
                            #     subpl2.plot(bin_mean, hist_to_plot_fin_ref, 'k--', label=)
                            # elif r == 'TRMM':
                            #     subpl2.plot(bin_mean, hist_to_plot_fin_ref, color='r', label=r)
                            # else:

                            subpl2.plot(bin_mean, hist_to_plot_fin_ref,
                                        color=colorlwdict[labref][0], linewidth=1,
                                        linestyle=colorlwdict[labref][2], label=labref)
                            # subpl2.plot(bin_mean, hist_to_plot_fin_ref, color='k', label=r[:-7])
                            if nb_freq_per_day == 1:
                                ylim_top = 0.41
                                ylim_bot = -0.05
                                subpl2.set_ylim([ylim_bot, ylim_top]) # -0.025, 0.025
                            # elif nb_freq_per_day == 8:
                            #     # WAfrica 3h:
                            #     subpl2.set_ylim([-0.009, 0.009])
                            else:
                                # Europe 1h:
                                # subpl2.set_ylim([-0.0015, 0.0025])
                                # W Africa 3h:
                                ylim_bot = -0.0035
                                ylim_top = 0.0061
                                subpl2.set_ylim([ylim_bot, ylim_top])
                            if (r != diff_res_todo[nb_r]):
                            #if r in ['ESM-r1',]:
                            #if r not in ['TRMM', 'CMORPH', 'CMORPH1', 'GPCC', 'GPCP']:
                                shift_res = {'N512r': -0.0006, 'N512r_future': -0.0006, 'CP4A': 0.0001,
                                             'CP4A_future': 0.0001, }
                                shift = 0
                                if r in shift_res:
                                    shift = shift_res[r]
                                plot_percent_change_between_centiles(subpl2, hist_to_plot_fin_ref, hist_to_plot_fin,
                                                                     bin_mean,
                                                                     color=colorlwdict[r][0], between=(
                                                                     one_third, two_thirds),
                                                                     ylim_max=x_lim[1],
                                                                     shift=shift, ylim_top=ylim_top)
                                # print(' p99.99 statistics for {}'.format(r))
                                # plot_percent_change_between_centiles(subpl2, hist_to_plot_fin_ref, hist_to_plot_fin,
                                #                                      bin_mean,
                                #                                      color=colorlwdict[r][0], between=(
                                #         centiles_ref['p99'], centiles_ref['p99.99']),
                                #         # one_third, two_thirds),
                                #                                      ylim_max=x_lim[1],
                                #                                      shift=shift, plot_only='upper',
                                #                                      ylim_top=ylim_top)
                                # plot_fraction_and_break_pts(subpl2, hist_to_plot_fin_ref, hist_to_plot_fin, bin_mean,
                                #                             color=colorlwdict[r][0])

                            # if nb_freq_per_day ==1:
                            #    subpl2.set_ylim([-0.008, 0.008])
                            # else:
                            #    subpl2.set_ylim([-0.05, 0.05])
                            # centils = np.array(np.sort(centiles.values()))
                            # print 'centils :', centils
                            # centils_ref = np.array(np.sort(centiles_ref.values()))
                            # if r in '2p2f_on_n512':
                            #     centils = centils * 3600
                            #     centils_ref = centils_ref * 3600
                            # print 'centils :', centils_ref
                            shift = ylim_top * 0.1  # 0.0004 #for Europe
                            if (centiles_ref['p99.99'] - centiles['p99.99']) == 0:
                                pass
                                # subpl2.plot(centiles_ref['p99.99'], 0.7*ylim_top, 's',
                                #     color=colorlwdict[r][0],
                                #     linewidth=colorlwdict[r][1],
                                #     markersize=5,)
                                # subpl2.text(centiles_ref['p99.99']+x_lim[1]*0.02, 0.71*ylim_top,
                                #             '$p99.99^{th}$', fontsize=7)
                            # elif r not in ['TRMM', 'CMORPH', 'CMORPH1']:
                            #     centile_yplace_on_graph = 0.75*ylim_top # 0.004
                            #     subpl2.arrow(centiles_ref['p99.999'], centile_yplace_on_graph - shift * nb_r,  # 0.0005 for Europe
                            #                  centiles['p99.999'] - centiles_ref['p99.999'], 0,
                            #                  width=0.00005,
                            #                  head_width=0.0001, head_length=1,
                            #                  fc=colorlwdict[r][0], ec=colorlwdict[r][0])
                            #
                            #     value = (centiles['p99.999'] - centiles_ref['p99.999']) / centiles_ref['p99.999'] * 100.
                            #     sign = value / np.abs(value)
                            #     add_sign = {-1: '', 1: '+'}
                            #     centile_text_yplace_on_graph = 0.75*ylim_top + ylim_top * 1/60.
                            #     subpl2.text(centiles_ref['p99.999'], centile_text_yplace_on_graph - shift * nb_r,  # 0.0007 for Europe
                            #                 '{}{:.0f}%'.format(add_sign[sign], value),
                            #                 fontsize=7, color=colorlwdict[r][0])
                            #     centile_yplace_on_graph = 0.45*ylim_top # 0.0028
                            #     subpl2.arrow(centiles_ref['p99.99'], centile_yplace_on_graph - shift * nb_r,  # 0.0005 for Europe
                            #                  centiles['p99.99'] - centiles_ref['p99.99'], 0,
                            #                  width=0.00005,
                            #                  head_width=0.0001, head_length=1,
                            #                  fc=colorlwdict[r][0], ec=colorlwdict[r][0])
                            #     value = (centiles['p99.99'] - centiles_ref['p99.99']) / centiles_ref[
                            #                         'p99.99'] * 100.
                            #     sign = value / np.abs(value)
                            #     add_sign = {-1: '', 1: '+'}
                            #     centile_text_yplace_on_graph = centile_yplace_on_graph + ylim_top * 1/60.
                            #     subpl2.text(centiles_ref['p99.99'], 0.0029 - shift * nb_r,  # 0.0007 for Europe
                            #                 '{}{:.0f}%'.format(add_sign[sign], value),
                            #                 fontsize=7, color=colorlwdict[r][0])

                    elif plot_type == 'centiles':
                        centils = np.array(np.sort(centiles.values()))
                        shift_dict = {'N512r': 0, 'CP4A': 0.8, 'N512r_future': 0, 'CP4A_future': 0.8, 'CMORPH': 0.4,
                                      'TRMM': 0.4}
                        abscissis = np.arange(0 + shift_dict[r], len(centils) * 1.5, 1.5)
                        # if r is not 'CMORPH':
                        print 'correcting mm/h into mm/3h: x3'
                        print r, centils * 3
                        subpl2.plot(abscissis, centils * 3, 's',
                                    color=colorlwdict[r][0],
                                    linewidth=colorlwdict[r][1],
                                    markersize=8,
                                    label='{}'.format(r))
                        # else:
                        #    print r, centils
                        #    subpl2.plot(abscissis, centils, 's',
                        #            color = colorlwdict[r][0], 
                        #            linewidth = colorlwdict[r][1], 
                        #            markersize = 8,
                        #            label='{}'.format(r))
                        arrow_dict[r] = centils
                        textt = ['{:.1f}'.format(n) for n in np.sort(centiles.values())]
                        # print textt
                        # plt.text([0.2,1.2,2.2,3.2], np.sort(centiles.values()),textt,
                        #         color = colorlwdict[r][0])
                        ylabel = 'precipitation rate (mm/{}h)'.format(int(24. / nb_freq_per_day))
                        subpl2.set_xticks(np.arange(0 + 0.4, len(centils) * 1.5, 1.5))
                        subpl2.set_xticklabels(
                            [r'$99^{th}$', r'$99.5^{th}$', r'$99.9^{th}$', r'$99.99^{th}$', r'$99.999^{th}$'],
                            rotation=40)
                        plt.xlim((-0.5, 7.3))
                    else:
                        subpl2.plot(bin_mean, hist_to_plot_fin, '-',
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1], linestyle=colorlwdict[r][2],
                                    label='{}'.format(r))
                        ylabel = 'Weighted frequency'
                        subpl2.set_ylim([0.000001, 1])
            # if plot_type == 'frac':
            if (plot_type == 'centiles') & ('N512r' in arrow_dict.keys()):
                for num, xabs in enumerate(np.arange(0, len(centils) * 1.5, 1.5)):
                    subpl2.arrow(xabs + shift_dict['N512r'], arrow_dict['N512r'][num] * 3,
                                 0, arrow_dict['N512r_future'][num] * 3 - arrow_dict['N512r'][num] * 3,
                                 head_width=0.2, head_length=0.3,
                                 fc=colorlwdict['N512r_future'][0], ec=colorlwdict['N512r_future'][0])
                    subpl2.arrow(xabs + shift_dict['CP4A'], arrow_dict['CP4A'][num] * 3,
                                 0, arrow_dict['CP4A_future'][num] * 3 - arrow_dict['CP4A'][num] * 3,
                                 head_width=0.2, head_length=0.3,
                                 fc=colorlwdict['CP4A_future'][0], ec=colorlwdict['CP4A_future'][0])
                    subpl2.text(xabs + shift_dict['N512r'],
                                (arrow_dict['N512r_future'][num] * 3 - arrow_dict['N512r'][num] * 3) / 2. +
                                arrow_dict['N512r'][num] * 3 - 1.0,
                                '+{:.0f}\%'.format(
                                    (arrow_dict['N512r_future'][num] - arrow_dict['N512r'][num]) / arrow_dict['N512r'][
                                        num] * 100.),
                                fontsize=7, color=colorlwdict['N512r_future'][0])
                    subpl2.text(xabs + shift_dict['CP4A'],
                                (arrow_dict['CP4A_future'][num] * 3 - arrow_dict['CP4A'][num] * 3) / 2. +
                                arrow_dict['CP4A'][num] * 3 - 1.0,
                                '+{:.0f}\%'.format(
                                    (arrow_dict['CP4A_future'][num] - arrow_dict['CP4A'][num]) / arrow_dict['CP4A'][
                                        num] * 100.),
                                fontsize=7, color=colorlwdict['CP4A_future'][0])
                    subpl2.bar(np.arange(0, len(centils) * 1.5, 1.5) + 0.4, np.ones(len(centils)) * 180 * 1.10,
                               # (np.max(arrow_dict['N512r_future'])*3*1.10),
                               np.ones(len(centils)) * 1.2, alpha=0.05, color='grey', align='center', ec=None)
                    plt.ylim([0, 180])  # np.max(arrow_dict.values()[-1])*3*1.05])
            if (np.mod(numplot - 1, numcol) == 0):
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='on')
                subpl2.set_ylabel(ylabel)
            elif (np.mod(numplot, numcol) == 0):
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='off',
                    left='off',
                    labelright='on',
                    right='on')
                subpl2.set_ylabel(ylabel)
            else:
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='off')  # 'off')
                subpl2.set_ylabel(ylabel)
            if not plot_type in ['frac', 'contrib', 'centiles']:
                if not diff_btw_res:
                    subpl2.set_yscale('log')
            if (numplot) > (num_lines - 1) * numcol:
                # if not iris.util.approx_equal(self.bins[2]-self.bins[1], self.bins[21]-self.bins[20]):
                #   subpl2.set_xlabel('PPN (mm/day)')
                # else:
                if plot_type == 'centiles':
                    subpl2.set_xlabel('centiles')
                    subpl2.set_ylabel('(mm/{}h)'.format(int(24. / nb_freq_per_day)))
                else:
                    subpl2.set_xlabel('bin rate (mm/{}h)'.format(int(24. / nb_freq_per_day)))
                    # subpl2.set_xlim([0, x_lim])
            if not plot_type == 'centiles':
                subpl2.set_xscale('log')
                subpl2.set_xlim([x_lim[0], x_lim[1]])
            # subpl2.grid(True)
            if 'json' in reg:
                reg = reg.split('_')[1:-1]
            if (diff_btw_res is not None) & (plot_type == 'frac'):
                subpl2.set_title(
                    ' '.join(w for w in reg))  # +', '+', '. join('{:2.2f}'.format(b) for b in bin_oc_tot))
            else:
                print bigregion, num_plots
                if (bigregion == (-20, 12, 0, 25)) & (num_plots == 5):
                    dict_nb = {'SEN': 'b. ', 'NSA': 'c. ', 'GSL': 'd. ', 'CSA': 'e. ', 'GOG': 'f.'}
                    dict_land = {'SEN': ', land only', 'NSA': '', 'GSL': ', land only', 'CSA': '', 'GOG': ', land only'}
                    if reg == 'CSA':
                        subpl2.set_title('{}{}{}'.format(dict_nb[reg], 'SSZ', dict_land[reg]), weight='bold')
                    else:
                        subpl2.set_title('{}{}{}'.format(dict_nb[reg], reg, dict_land[reg]), weight='bold')
                elif (bigregion == (-20, 12, 0, 25)) & (num_plots == 3):
                    dict_nb = {'SEN': 'c. ', 'NSA': 'b. ', 'GSL': 'd. ', 'CSA': 'd. ', 'GOG': 'f.'}
                    dict_land = {'SEN': ', land only', 'NSA': '', 'GSL': ', land only', 'CSA': '', 'GOG': ', land only'}
                    subpl2.set_title('{}{}{}'.format(dict_nb[reg], reg, dict_land[reg]), weight='bold')
                else:
                    if col_order is not None:
                        subpl2.set_title(' '.join(w for w in reg), color=col_order['_'.join(w for w in reg)])
                    else:
                        subpl2.set_title(' '.join(w for w in reg))
            if subregion_files:
                if col_order is not None:
                    print 'col_order', col_order
                    print 'for reg', reg, ' ', col_order['_'.join(w for w in reg)]
                    for ax_pos in ['bottom', 'top', 'right', 'left']:
                        subpl2.spines[ax_pos].set_color(col_order['_'.join(w for w in reg)])
                        subpl2.spines[ax_pos].set_linewidth(2)
            if regcount == 1:
                 subpl2.legend(prop={'size': 10}, loc='best')
            # legend_axes = plt.gcf().add_axes([0.56, 0.56, 0.2, 0.2])
            # subpl2.legend(prop={'size':9},bbox_to_anchor=(-0.6,0.45,0.2,0.55), loc=4, borderaxespad=0.1)
            # for group in pie:
            #   for x in group:
            #       x.set_visible(False)
            # subpl2.legend(loc='upper left', prop={'size':10})
            # plt.legend(subpl2, legend_axes, prop={'size':8})
        plt.show()
        plt.savefig(target + '.png', dpi=300)
        plt.savefig(target + '.pdf')

####update11/10
    def save_subplot_with_spread(self, reg_todo, plot_type, x_lim, colorlwdict, figdir, average=False, season=None, bigregion=None,
                     subregion_files=None, nb_freq_per_day=1., dict_names={}, dict_ensemble=None,
                     bootstrap_config=None, bootstrap_or_centiles='bootstrap', percentiles=[2.5, 97.5], sig_diff_btw_list=None, sig_val=0.01, vals=[], plot=True):
        '''
        Save all the plots

        :param dict reg_todo: Contains all the region shapes
        :param str plot_type: Type of data being plotted to be included in filename: contrib, frac, normed, normed_wet, centiles
        :param int x_lim: Sets the max of the x-axis
        :param dict colorlwdict: contains [ line color, linewidth, linetype] for each model
        :param str figdir: directory where the figure will be saved
        :param bool average: not used anymore, used to be to take the average of different ensemble members
        :param str season: e.g 'jja'
        :param dict subregion_files: for interface_plot_european_masked_subregions: uses masks to define regions
        :param float nb_freq_per_day: 24. if hourly, 8. if 3 hourly, 1. if daily
        :param list diff_btw_res: diff_btw_res = ['dataset_ref', 'dataset1', 'dataset_ref', 'dataset2',]
                                  will plot dataset1 - dataset_ref and dataset2 - dataset_ref
        :param dict dict_names: optional, if you want a different label name than the resolution name
        :param dict dict_ensemble: e.g. {'PRIMAVERA': [model1, model2, model3, ...], } will plot mean + 25th/75th centiles
        :param interannual_spread: will plot mean + inter-annual spread around mean.
        '''
        from matplotlib import style
        style.use('seaborn-colorblind')
        freqs = {'1': 'd', '8': '3h', '24': 'h', '96': '15mn'}
        res_todo = self.resolutions
        hfiter = False
        season_print = season
        if isinstance(season, list):
            season_print = ['_'.join(s for s in season)]
        target = "{}/{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season_print,
                                             int(24. / nb_freq_per_day))
        regions = reg_todo
        if reg_todo == 'masked':
            hfiter = True
            print self.hists.keys()
            reg_todo = self.hists[self.hists.keys()[0]]['masked'].keys()
            regions = None
            target = "{}/{}_{}_{}_{}_{}_{}".format(figdir, '_'.join(s[4:-5] for s in reg_todo), res_todo[0], plot_type,
                                                   nb_freq_per_day, freqs[str(int(nb_freq_per_day))], season_print)
      
        print target
        if plot:
         num_plots = len(reg_todo)
         # ### plots the regions of interest first ### does not plot the regions if only 1 region.
         fig, subpl, numcol, start_plot, num_lines, col_order = set_up_fig(num_plots, reg_todo=regions,
                                                                          bigregion=bigregion,
                                                                          subregion_files=subregion_files,
                                                                          season=season)

         if regions:
            if isinstance(bigregion, tuple):
                if subpl:
                    subpl.set_extent(bigregion, ccrs.PlateCarree())  # Africa
        percent_sig = {}
        for regcount, reg in enumerate(reg_todo):
            print 'reg', reg
            reg_short = reg
            if len(reg.split('_')) > 1:
                reg_short = reg.split('_')[1]
            percent_sig[reg_short] = {}
            if plot:
             numplot = start_plot + regcount
            # new subplot for each region
             subpl2 = fig.add_subplot(num_lines, numcol, numplot)
            if season:
                season = reg
            print reg
            if reg is not None:
                bin_oc_tot = []
                arrow_dict = {}
                dict_all_hist_to_plot = {}
                for nb_r, r in enumerate(res_todo):
                    if iris.util.approx_equal(self.bins[2] - self.bins[1],
                                              self.bins[21] - self.bins[20]):  # traditional distrib with regular bins
                        bins1 = self.bins * 86400. / nb_freq_per_day
                    else:  # Martin distribution
                        #print 'Martin distribution, assuming given in kg/m2/s, check in process_data.py', self.bins.shape
                        bins1 = self.bins * 86400. / nb_freq_per_day
                    # if nb_r==0:
                    #     for xset in bins1:
                    #         plt.axvline(x=xset, linewidth=0.3, color='k')
                    if average:
                        #		runs=self.config['histfile'][r]['average']
                        hf = 'average'
                    else:
                        bin_mean = (bins1[:-1] + bins1[1:]) / 2.
                        #print 'bin mean', bin_mean
                        ####### get the histograms to plot #########
                        hist_to_plot, _, _, _ = self.get_hist2plot(r, reg, hfiter)
                        interannual_spread = True
                        for ens in dict_ensemble.keys():
                            if r in dict_ensemble[ens]:
                                interannual_spread = False
                        if interannual_spread:
                            print 'interannual spread to be calculated'
                            hist_to_plot_dict = hist_to_plot[0]
                        else:
                            print 'intermember spread to be calculated'
                            hist_to_plot_dict = {}
                            hist_to_plot_dict['0'] = get_mean_10_year_change(hist_to_plot[0])
                        dict_all_hist_to_plot[r] = {}
                        for year in hist_to_plot_dict.keys():
                            hist_to_plot = hist_to_plot_dict[year]
                            if plot_type == 'normed':
                                hist_to_plot_fin = hist_to_plot / hist_to_plot.sum()
                            elif plot_type == 'normed_wet':
                                if nb_freq_per_day == 1.:
                                    hist_to_plot_fin = hist_to_plot[1:] / hist_to_plot[1:].sum()
                                    bins_to_plot = bin_mean[1:]
                                else:
                                    hist_to_plot_fin = hist_to_plot[2:] / hist_to_plot[2:].sum()
                                    bins_to_plot = bin_mean[2:]
                            elif plot_type == 'frac':
                                hist_to_plot_fin = hist_to_plot * bin_mean / sum(hist_to_plot * bin_mean) * 100
                            elif plot_type == 'contrib':
                                hist_to_plot_fin = hist_to_plot * bin_mean / hist_to_plot.sum()
                                #print 'hirst_to_plot_fin', hist_to_plot_fin
                            elif plot_type == 'centiles':
                                centiles = self.get_centiles(r, reg, hfiter)
                            else:
                                hist_to_plot_fin = hist_to_plot
                            dict_all_hist_to_plot[r][year] = hist_to_plot_fin

                if dict_ensemble is not None:
                    dict_all_hist_to_plot_new = {}
                    list_in_ens = []
                    for ens_name in dict_ensemble:
                        print ens_name
                        dict_all_hist_to_plot_new[ens_name] = {}
                        for num, r in enumerate(dict_ensemble[ens_name]):
                            list_in_ens.append(r)
                            ### give a fake year to ensemble members
                            print bootstrap_config[ens_name]
                            dict_all_hist_to_plot_new[ens_name][str(num+bootstrap_config[ens_name]['master_ref_year'])] = dict_all_hist_to_plot[r]['0']

                    print list_in_ens
                    for r in res_todo:
                        if r not in list_in_ens:
                            dict_all_hist_to_plot_new[r] = dict_all_hist_to_plot[r]
                    dict_all_hist_to_plot_fin = dict_all_hist_to_plot_new
                else:
                    dict_all_hist_to_plot_fin = dict_all_hist_to_plot


                ########### now plot, you may want to change some things in here ##############

                for r in dict_all_hist_to_plot_fin.keys():
                    if 'obs' in r:
                        if 'scale' in r:
                           bootstrap_config_to_do = bootstrap_config['obs_cordex50'][reg]
                        else:
                           bootstrap_config_to_do = bootstrap_config[r][reg]
                    else:
                        bootstrap_config_to_do = bootstrap_config[r]
                    print 'percentiles', percentiles
                    low_high_bounds, mean_to_plot = get_histo_significance_interval(dict_all_hist_to_plot_fin[r].copy(), r,
                                                                                        bootstrap_config_to_do, bin_mean,
                                                                                        percentiles, bootstrap_or_centiles)
                    #print 'mean_to_plot.data', mean_to_plot.data
                    if plot:
                     print 'ens. length', r, ' ', len(dict_all_hist_to_plot_fin[r])
                     ### plot individual curves
                     #for his in dict_all_hist_to_plot_fin[r]:
                     #   subpl2.plot(bin_mean, dict_all_hist_to_plot_fin[r][his],
                     #               color=colorlwdict[r][0], linewidth=0.5, linestyle=colorlwdict[r][2],
                     #               )

                     if bootstrap_or_centiles == 'bootstrap':
                       subpl2.plot(bin_mean, mean_to_plot.data,
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1], linestyle=colorlwdict[r][2],
                                    label='{}'.format(r))
                       low_bound = low_high_bounds.extract(iris.Constraint(percentile_over_bootstrap_sample_number=lambda p: p == percentiles[0]))
                       high_bound = low_high_bounds.extract(iris.Constraint(percentile_over_bootstrap_sample_number=lambda p: p == percentiles[1]))
                     else:
                       print low_high_bounds.coord('percentile_over_years') 
                       median = low_high_bounds.extract(iris.Constraint(percentile_over_years=lambda p: p == 50))
                       low_bound = low_high_bounds.extract(iris.Constraint(percentile_over_years=lambda p: p == 25))#percentiles[0]))
                       high_bound = low_high_bounds.extract(iris.Constraint(percentile_over_years=lambda p: p == 75))#percentiles[1]))
                       subpl2.plot(bin_mean, median.data,
                                    color=colorlwdict[r][0], linewidth=colorlwdict[r][1], linestyle=colorlwdict[r][2],
                                    label='{}'.format(r))
                     if not 'scale' in r:
                       subpl2.fill_between(bin_mean, low_bound.data, high_bound.data, color=colorlwdict[r][0], alpha=0.5)
                    
                    #low_bound = low_high_bounds.extract(iris.Constraint(percentile_over_bootstrap_sample_number=lambda p: p == 0.5))
                    #high_bound = low_high_bounds.extract(iris.Constraint(percentile_over_bootstrap_sample_number=lambda p: p == 99.5))
                   # subpl2.fill_between(bin_mean, low_bound.data, high_bound.data, color=colorlwdict[r][0], alpha=0.5)
                    
                    if plot:
                     if plot_type == 'normed':
                        print 'total sum of frequencies for {}, {}: {}'.format(r, reg, hist_to_plot_fin.sum())
                        ylabel = 'Normed frequency'
                        subpl2.set_ylim([0.000001, 0.2])  # daily Africa
                        subpl2.set_ylim([0.0000001, 0.1])  # daily Africa
                        # subpl2.set_ylim([0.00000002, 0.1])#[0.0000001,0.005])#[0.001, 0.5])#[0.0000001,0.005])#[0.0000001, 1])
                     elif plot_type == 'normed_wet':
                        print 'total sum of frequencies for {}, {}: {}'.format(r, reg, hist_to_plot_fin.sum())
                        ylabel = 'Normed frequency (wet hours)'
                        subpl2.set_ylim([0.000001, 1])
                     elif plot_type == 'frac':
                        ylabel = 'Fractional contribution (\%) '
                        subpl2.set_ylim([0, 3])
                            # plt.gca().set_yticks(range(4))
                            # plt.gca().set_yticklabels(range(4))
                     elif plot_type == 'contrib':
                        # if r is not 'TRMM':
                        lab = r
                        if r == 'N512r':
                            lab = 'R25'
                        if lab in dict_names.keys():
                            lab = dict_names[lab]
                        # if r is not 'TRMM':
                        if r in ['SEA', 'AMZ']:
                            subpl2.set_ylim([0, 0.21])
                        else:
                            subpl2.set_ylim([0, 0.23])  # 0.012]) #0.17]) #0.023])
                        if 24. / nb_freq_per_day == 1:
                            ylabel = 'Contribution (mm/h)'
                        elif nb_freq_per_day == 1:
                            ylabel = 'Contribution (mm/day)'
                        else:
                            ylabel = 'Contribution (mm/{}h)'.format(int(24. / nb_freq_per_day))
                     else:
                        ylabel = 'Weighted frequency'
                        subpl2.set_ylim([0.000001, 1])

                if sig_diff_btw_list:
                    aa = 0
                    for aa, sig_diff_btw in enumerate(sig_diff_btw_list):
                        name_diff = '{}vs{}'.format(sig_diff_btw[0], sig_diff_btw[1])
                        percent_sig[reg_short][name_diff] = {}
                        bootstrap_conf = {}
                        for r in sig_diff_btw:
                            if 'obs' in r:
                                bootstrap_conf[r] = bootstrap_config[r][reg]
                            else:
                                bootstrap_conf[r] = bootstrap_config[r]
                        p_val = test_hist_change_significance(dict_all_hist_to_plot_fin[sig_diff_btw[0]],
                                                              dict_all_hist_to_plot_fin[sig_diff_btw[1]],
                                                              sig_diff_btw[0],
                                                              sig_diff_btw[1],
                                                              bootstrap_conf[sig_diff_btw[0]],
                                                              bootstrap_conf[sig_diff_btw[1]], bin_mean
                                                              )
                        if plot:                                      
                         if sig_diff_btw == 1: #['CORDEX-44', 'PRIMAVERA', ]: #only plot diff between first two datasets
                           ax2 = subpl2.twinx()
                           ax2.plot(bin_mean, p_val, '+',
                                    color='k', alpha=0.4)
                           ax2.set_ylim([0, 0.5])
                           ax2.axhline(0.01, linewidth=1)
                           fig.canvas.draw()
                           ax2.yaxis.set_label_position("right")
                           ax2.tick_params(axis="y",direction="in", pad=-24, color='grey')
                           labels = [item.get_text() for item in ax2.get_yticklabels()]
                           labels[0] = ''
                           labels[-1] = '' ; labels[-2] = ''
                           ax2.set_yticklabels(labels, color='grey')

                        for val_min_max in vals:
                            p_val_masked = np.ma.masked_where(bin_mean < val_min_max[0], p_val.copy())
                            p_val_masked = np.ma.masked_where(bin_mean > val_min_max[1], p_val_masked)
                            num_vals = p_val_masked.count()
                            p_val_masked_sig = np.ma.masked_where(p_val_masked > sig_val, p_val_masked)
                            percent_sig[reg_short][name_diff]['{}to{}'.format(val_min_max[0], val_min_max[1])] = p_val_masked_sig.count() / np.float(num_vals)

            if plot:
             if (np.mod(numplot - 1, numcol) == 0):
                subpl2.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='on')
                #if sig_diff_btw_list:
                #    ax2.tick_params(axis='y', labelright='off', right='off')
                subpl2.set_ylabel(ylabel)
             elif (np.mod(numplot, numcol) == 0):
                subpl2.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='off',
                    left='off',
                    labelright='on',
                    right='on')
                subpl2.set_ylabel(ylabel)
             else:
                subpl2.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='off')  # 'off')
                #if sig_diff_btw_list:
                #    ax2.tick_params(axis='y', labelright='off', right='off')
                subpl2.set_ylabel(ylabel)
             if (numplot) > (num_lines - 1) * numcol:

                if plot_type == 'centiles':
                    subpl2.set_xlabel('centiles')
                    subpl2.set_ylabel('(mm/{}h)'.format(int(24. / nb_freq_per_day)))
                else:
                    subpl2.set_xlabel('bin rate (mm/{}h)'.format(int(24. / nb_freq_per_day)))
                    # subpl2.set_xlim([0, x_lim])
             if not plot_type == 'centiles':
                subpl2.set_xscale('log')
                subpl2.set_xlim([x_lim[0], x_lim[1]])
             # subpl2.grid(True)
             if 'json' in reg:
                reg = reg.split('_')[1:-1]
             print bigregion, num_plots
             if (bigregion == (-20, 12, 0, 25)) & (num_plots == 5):
                dict_nb = {'SEN': 'b. ', 'NSA': 'c. ', 'GSL': 'd. ', 'CSA': 'e. ', 'GOG': 'f.'}
                dict_land = {'SEN': ', land only', 'NSA': '', 'GSL': ', land only', 'CSA': '', 'GOG': ', land only'}
                if reg == 'CSA':
                    subpl2.set_title('{}{}{}'.format(dict_nb[reg], 'SSZ', dict_land[reg]), weight='bold')
                else:
                    subpl2.set_title('{}{}{}'.format(dict_nb[reg], reg, dict_land[reg]), weight='bold')
             elif (bigregion == (-20, 12, 0, 25)) & (num_plots == 3):
                dict_nb = {'SEN': 'c. ', 'NSA': 'b. ', 'GSL': 'd. ', 'CSA': 'd. ', 'GOG': 'f.'}
                dict_land = {'SEN': ', land only', 'NSA': '', 'GSL': ', land only', 'CSA': '', 'GOG': ', land only'}
                subpl2.set_title('{}{}{}'.format(dict_nb[reg], reg, dict_land[reg]), weight='bold')
             else:
                if col_order is not None:
                    subpl2.set_title(' '.join(w for w in reg), color=col_order['_'.join(w for w in reg)])
                else:
                    subpl2.set_title(' '.join(w for w in reg))
             if subregion_files:
                if col_order is not None:
                    print 'col_order', col_order
                    print 'for reg', reg, ' ', col_order['_'.join(w for w in reg)]
                    for ax_pos in ['bottom', 'top', 'right', 'left']:
                        subpl2.spines[ax_pos].set_color(col_order['_'.join(w for w in reg)])
                        subpl2.spines[ax_pos].set_linewidth(2)
             if regcount == 1:
                subpl2.legend(prop={'size': 10}, loc='best')
 
        if plot:
         plt.savefig(target + '.png', dpi=300)
         plt.savefig(target + '.pdf')

        return percent_sig


####end update11/10
    def get_hist2plot(self, r, reg, hfiter, d1_or_2d='1d'):
        if d1_or_2d == '2d':
            hist_list = [self.hists.copy(), self.hist2d.copy()]
        else:
            hist_list = [self.hists.copy(), ]
        for num, histo in enumerate(hist_list):
            if hfiter:
                print 'r', r, 'reg', reg
                hists_to_plot = histo[r]['masked'][reg]
                nb_years = self.end_year[r]['masked'][reg] - self.start_year[r]['masked'][reg] + 1
                nb_spatial_points = self.nb_spatial_points[r]['masked'][reg]
                totnum = self.totnum[r]['masked'][reg]
            else:
                hf = histo[r][reg].keys()[0]
                hists_to_plot = histo[r][reg][hf]
                nb_years = self.end_year[r][reg][hf] - self.start_year[r][reg][hf] + 1
                nb_spatial_points = self.nb_spatial_points[r][reg][hf]
                totnum = self.totnum[r][reg][hf]
            # print hists_to_plot.values()[0]
            if len(np.array(hists_to_plot.values()[0]).shape) == 1:
                for year in hists_to_plot:
                    if len(np.array(hists_to_plot[year])) == len(self.bins):
                        hists_to_plot[year] = np.array(hists_to_plot[year][0:-1])
                    else:
                        hists_to_plot[year] = np.array(hists_to_plot[year])
            elif len(np.array(hists_to_plot.values()[0]).shape) == 2:
                for year in hists_to_plot:
                    hists_to_plot[year] = np.array(hists_to_plot[year])
            else:
                for year in hists_to_plot:
                    hists_to_plot[year] = np.array(hists_to_plot[year])
            hist_list[num] = hists_to_plot
        return hist_list, nb_years, nb_spatial_points, totnum

    def get_yearly_time_series(self, r, reg, hfiter):
        if hfiter:
            print 'r', r, 'reg', reg
            time_series = self.yearly_time_series[r]['masked'][reg]
            year_start = self.start_year[r]['masked'][reg]
            year_end = self.end_year[r]['masked'][reg]
        else:
            hf = self.yearly_time_series[r][reg].keys()[0]
            time_series = self.yearly_time_series[r][reg][hf]
            year_start = self.start_year[r][reg][hf]
            year_end = self.end_year[r][reg][hf]
        year_range = np.arange(year_start, year_end + 1)
        return year_range, time_series

    def get_centiles(self, r, reg, hfiter):
        if hfiter:
            print 'r', r, 'reg', reg
            centiles = self.centiles[r]['masked'][reg]
        else:
            hf = self.centiles[r][reg].keys()[0]
            centiles = self.centiles[r][reg][hf]
        return centiles

    def save_time_series(self, reg_todo, plot_type, x_lim, colorlwdict, figdir, average=False, season=None,
                         bigregion=None,
                         subregion_files=None, nb_freq_per_day=1., diff_btw_res=None, dict_names={}):
        '''
        Save all the plots

        :param dict reg_todo: Contains all the region shapes
        :param int x_lim: Sets the max of the x-axis
        :param str plot_type: Type of data being plotted to be included in filename
        :param dict hist: Contains 1D numpy arrays being the histograms
        :param int num_data: Number of data sets being used
        :param bool normed: Whether the plot should be normalised or not
        :param float nb_freq_per_day: 24. if hourly, 8. if 3 hourly, 1. if daily
        :param dict dict_names: optional, if you want a different label name than the resolution name
        '''

        freqs = {'1': 'd', '8': '3h', '24': 'h', '96': '15mn'}
        print 'diff_btw_res', diff_btw_res, bool(diff_btw_res)
        if diff_btw_res:
            # plot the difference: diff_btw_res[1]-diff_btw_res[0], diff_btw_res[3]-diff_btw_res[2], ...
            res_todo = diff_btw_res[1::2]  # take the first and then every 2 datasets
            diff_res_todo = diff_btw_res[::2]  # take reference datasets every 2 datasets.
        else:
            res_todo = self.resolutions
        hfiter = False
        if diff_btw_res:
            target = "{}/{}_{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season,
                                                int(24. / nb_freq_per_day), diff_btw_res[0])
        else:
            target = "{}/{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season,
                                             int(24. / nb_freq_per_day))
        regions = reg_todo
        if reg_todo == 'masked':
            hfiter = True
            reg_todo = self.hists[self.hists.keys()[0]]['masked'].keys()
            regions = None
            target = "{}/{}_{}_{}_{}_{}_{}".format(figdir, '_'.join(s[4:-5] for s in reg_todo), res_todo[0], plot_type,
                                                   nb_freq_per_day, freqs[str(int(nb_freq_per_day))], season)
        #   if (len(reg_todo) == 1) & (len(frequencies)>1):
        ################ adddddd ############
        print target
        num_plots = len(reg_todo)
        fig, subpl, numcol, start_plot, num_lines, col_order = set_up_fig(num_plots, reg_todo=regions,
                                                                          bigregion=bigregion,
                                                                          subregion_files=subregion_files)
        # if reg_todo:
        #     for r in srex_reg.keys():
        #         x, y = srex_reg[r].exterior.xy
        #         subpl.plot(x, y, color='k')  # plot regional box lines
        #         xc, yc = np.array(srex_reg[r].exterior.centroid)
        #         subpl.text(xc, yc, r, horizontalalignment='center', verticalalignment='center', color='k',
        #                    fontdict={'weight': 'bold'})  # add region codes to boxes

        if regions:
            if len(regions) > 1:
                subpl.set_extent(bigregion, ccrs.PlateCarree())  # Africa
        for regcount, reg in enumerate(reg_todo):
            numplot = start_plot + regcount
            subpl2 = fig.add_subplot(num_lines, numcol, numplot)
            if season:
                season = reg
            print reg
            if reg is not None:
                bin_oc_tot = []
                arrow_dict = {}
                for nb_r, r in enumerate(res_todo):
                    shift = np.float(nb_r) * 0.9 / len(res_todo)
                    if diff_btw_res:
                        year_range, time_series = self.get_yearly_time_series(r, reg, hfiter)
                        year_range_ref, ref_time_series = self.get_yearly_time_series(diff_res_todo[nb_r], reg, hfiter)
                        if year_range == year_range_ref:
                            subpl2.bar(year_range + shift, time_series - ref_time_series, color=colorlwdict[r][0],
                                       alpha=0.8)
                    else:
                        year_range, time_series = self.get_yearly_time_series(r, reg, hfiter)
                        if r == 'CP4A_future':
                            year_range = np.array((1997, 1998, 1999, 2000, 2002, 2003, 2004, 2005))
                        print r, len(year_range + shift), len(time_series)
                        subpl2.bar(year_range + shift, time_series, color=colorlwdict[r][0], alpha=0.8, label=r,
                                   width=0.2)

            ylabel = 'total precipitation'
            if (np.mod(numplot - 1, numcol) == 0):
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='on')
                subpl2.set_ylabel(ylabel)
            elif (np.mod(numplot, numcol) == 0):
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='off',
                    labelright='on')
                subpl2.set_ylabel(ylabel)
            else:
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='on')  # 'off')
                # subpl2.set_ylabel(ylabel)
            plt.title(reg)
            if regcount == 4:
                subpl2.legend(prop={'size': 6}, loc='best')
        plt.show()
        plt.savefig(target + '.png', dpi=300, bbox_inches='tight')

    def save_2dhists(self, wet_or_dry_spells, thrs, reg_todo, x_lim, ylim, season=None, bigregion=None,
                     subregion_files=None, hist=None, freq_factor=1., replace_names=False, diff_btw_res=None,
                     integral_over_intensity=False, colorlwdict=None, numref_list=[1, 3], no_map=False, resol={},
                     figdir='', bootstrap_config=None):
        '''
        Save all the plots
        :param dict reg_todo: Contains all the region shapes
        :param int x_lim: Sets the max of the x-axis
        :param dict hist: Contains 1D numpy arrays being the histograms
        :param int num_data: Number of data sets being used
        :param bool normed: Whether the plot should be normalised or not
        :param freq_factor: number of samples per day at that frequency (mm/h -> 24.)
        :param tuple bigregion: (lon_min, lon_max, lat_min, lat_max)
        :param numref the plot number where the reference dataset is plotted (must be consistent with diff_btw_res,
        :param no_map if numplot > 2 and you don't want to plot the region map.
        :return str/int/float/iris.cube/ndarray 
         '''
        dict_unit = {'24': 'h', '1': 'days'}
        # resol = {
        #     'germany_on_12km': 'OBS on 12km',
        #     'switzerland_on_12km': 'OBS on 12km',
        #     'nimrod_on_12km': 'OBS on 12km',
        #     '2p2erai_on_12km': 'UKMO 2.2km',
        #     'mi-ao438': 'UKMO 12km',
        #     'mi-ao438_38': 'UKMO 12km',
        #     '2p2eth_on_12km': 'ETH 2.2km',
        #     '2p2eth38_on_12km': 'ETH 2.2km',
        #     '12eth_on_12km': 'ETH 12km',
        #     '12eth38_on_12km': 'ETH 12km',
        #     # 'N512r':'d. R25',
        #     # 'CP4A':'b. CP4A',
        #     'CMORPH':'d. CMORPH',
        #     'TRMM':'a. TRMM',
        # }

        plot_order = {'mi-ao438': 1, '2p2erai_on_12km': 2, 'germany_on_12km': 3, 'switzerland_on_12km': 3,
                      'nimrod_on_12km': 3, '12eth_on_12km': 4, '2p2eth_on_12km': 5}

        res_todo = self.resolutions
        if integral_over_intensity:
            num_plots = 1
        else:
            num_plots = len(res_todo)
        print 'num_plots', num_plots
        diff_name = 'ref'
        if diff_btw_res:
            diff_name = 'diff'

        if not subregion_files: # region defined at the start of this file as a polygon + land-sea mask (land-points only)
            hfiter = False
            fig, subpl, numcol, start_plot, num_lines, _ = set_up_fig(num_plots, [reg_todo, ],
                                                                      no_map=no_map, numref_list=numref_list)
            if subpl:
                x, y = srex_reg[reg_todo].exterior.xy
                subpl.plot(x, y, color='k')  # plot regional box lines
                xc, yc = np.array(srex_reg[reg_todo].exterior.centroid)
                subpl.text(xc, yc, reg_todo, horizontalalignment='center', verticalalignment='center', color='m',
                           fontdict={'weight': 'bold'})
                subpl.set_extent(bigregion, ccrs.PlateCarree())
            target = "{}/2dhists_{}_{}_{}_{}_{}_{}.png".format(figdir, '_'.join(s for s in res_todo), reg_todo,
                                                            wet_or_dry_spells, thrs, season, diff_name)
        else: # masked regions as input (e.g. country) in subregion_files
            hfiter = True
            fig, subpl, numcol, start_plot, num_lines, col_order = set_up_fig(num_plots, None, bigregion=bigregion,
                                                                              subregion_files=subregion_files,
                                                                              no_map=no_map, numref_list=numref_list)
            reg_todo = hist
            print reg_todo
            target = "{}/2dhists_{}_{}_{}_{}_{}_{}.png".format(figdir, '_'.join(s for s in res_todo), reg_todo[4:-5],
                                                            wet_or_dry_spells, thrs, season, diff_name)
        res_ref = res_todo[0]
        if diff_btw_res:
            res_ref_list = diff_btw_res[::2] # reference list
            res_todo = diff_btw_res[1::2]
        for rescount, res in enumerate(res_todo):
            # if rescount == 0:
            #    numcol = numcol
            # else:
            # numcol = 3
            # if start_plot+rescount < 3:
            #    numplot=start_plot+rescount+
            #    subpl2 = fig.add_subplot(num_lines, numcol, numplot)
            # else:
            if replace_names:
                numplot = plot_order[res] + 1
            else:
                numplot = start_plot + rescount  # +4-numcol
            if not ((integral_over_intensity) and (rescount > 0)):
                print num_lines, numcol, numplot
                subpl2 = fig.add_subplot(num_lines, numcol, numplot)
            if not season:
                season = reg_todo
            print res

            if np.mean(self.bins) > 20:
                bins1 = self.bins
            else:
                bins1 = self.bins * 86400. / freq_factor
                # for hf in runs:

            plotting_dur_bins = self.duration_bins[:-1] * 3 # *3for Africa 3h
            plotting_intens_bins = self.bins[:-1] * 86400. / freq_factor
            X, Y = np.meshgrid(plotting_dur_bins, plotting_intens_bins)

            if (diff_btw_res is not None) & (numplot not in numref_list):
                cmap = cm.get_cmap('BrBG', 11)
            else:
                cmap = cm.get_cmap('terrain_r', 12) #'brewer_RdYlGn_11')
                # cmap = cm.get_cmap('brewer_GnBu_09')

            hist, nb_years, nb_spatial_pts, totnum = self.get_hist2plot(res, reg_todo, hfiter, d1_or_2d='2d')
            hist = hist[0] # 2d hist
            # hist = np.transpose(hist)
            hists_to_plot = get_mean_10_year_change(hist.copy())
            totnum = get_mean_10_year_change(totnum.copy())
            if (res in ['N512r_future', 'CP4A', 'N512r', 'CP4A_future']) and (nb_years == 11):
                print 'number_years for {} changed to 10'.format(res)
                nb_years = 10.
            totnum_per_cell_per_year = np.float(totnum) / nb_years / nb_spatial_pts

            if diff_btw_res:
                print 'diff:{}-{}'.format(res, res_ref)
                res_ref = res_ref_list[rescount]
                hist_ref, nb_years_ref, nb_spatial_pts_ref, totnum_ref = self.get_hist2plot(res_ref, reg_todo, hfiter, d1_or_2d='2d')
                hist_ref = hist_ref[0] # 2d hist
                hists_to_plot_ref = get_mean_10_year_change(hist_ref.copy())
                totnum_ref = get_mean_10_year_change(totnum_ref.copy())
                if (res_ref in ['CP4A', 'N512r']) and (nb_years_ref == 11):
                    print 'number_years_ref for {} changed to 10'.format(res_ref)
                    nb_years_ref = 10.
                totnum_per_cell_per_year_ref = np.float(totnum_ref) / nb_years_ref / nb_spatial_pts_ref
                print 'nb years: {}, nb_years_ref: {}'.format(nb_years, nb_years_ref)
                print 'nb spatial points: {}, nb_spatial_points_ref: {}'.format(nb_spatial_pts, nb_spatial_pts_ref)
                print 'totnum:{}, totnum_ref:{}'.format(totnum, totnum_ref)

            if integral_over_intensity:
                hists_to_plot = np.sum(hists_to_plot, 0)
                if diff_btw_res:
                    hists_to_plot_ref = np.sum(hists_to_plot_ref, 0)
                plotting_intens_bins = None

            sig_mask = None
            if diff_btw_res:
                if (numplot not in numref_list):   # 4 for Europe
                    if bootstrap_config:
                        sig_mask = generate_sig_mask(hist_ref.copy(), hist.copy(), res_ref, res, bootstrap_config, plotting_intens_bins,
                                                 plotting_dur_bins)
                    hists_to_plot1 = (hists_to_plot - hists_to_plot_ref) / np.float(nb_spatial_pts_ref*nb_years_ref) # hists_to_plot/totnum - hists_to_plot_ref/totnum_ref #nb_years * nb_years_ref
                    hists_to_plot = hists_to_plot1
                    totnum_diff = (totnum_per_cell_per_year - totnum_per_cell_per_year_ref) / totnum_per_cell_per_year_ref * 100
                    # cube2plot = hist2cube(hists_to_plot, plotting_intens_bins, plotting_dur_bins)

            if integral_over_intensity:
                categories = ['$<$5days', '5-10days', '10-15 days', '$>$15 days']
                grouped_hist = [np.sum(hists_to_plot[0:5]), np.sum(hists_to_plot[5:10]), np.sum(hists_to_plot[10:15]),
                                np.sum(hists_to_plot[15:])]
                if res in resol.keys():
                    res_name = resol[res]
                else:
                    res_name = res
                cms = subpl2.bar(np.array([0, 1.2, 2.4, 3.6]) + 0.2 * (rescount + 1), grouped_hist, 0.2,
                                 color=colorlwdict[res][0], label=res_name)
                subpl2.set_xticks([0.6, 1.8, 3.0, 4.2])
                subpl2.set_xticklabels(categories, rotation=0, fontsize=20)
                # cms = subpl2.plot(self.duration_bins[:-1], hists_to_plot, label = res,
                #                  color = colorlwdict[res][0], linewidth = colorlwdict[res][1],)
            else:
                if (diff_btw_res is not None) and (numplot not in numref_list):  # 4 for Europe
                    if sig_mask is not None:
                        new_hists_to_plot = np.ma.array(hists_to_plot, mask=(sig_mask == False))
                        hists_to_plot = new_hists_to_plot
                    cmap.set_bad('white', 1.0)
                    cms = plt.pcolormesh(np.transpose(hists_to_plot),  # levels=[10**(-x) for x in range(1, 10)],
                                            cmap=cmap, norm=col.SymLogNorm(linthresh=0.001, linscale=0.1,
                                            vmin=-10, vmax=10))
                else:
                    cms = plt.pcolormesh(np.transpose(hists_to_plot) / np.float(nb_years*nb_spatial_pts),  # levels=[10**(-x) for x in range(1, 10)],
                                            cmap=cmap, norm=col.LogNorm(), vmin=10 ** (-2), vmax=10)
            if reg_todo == 'GOG':
                subpl2.set_xlabel('duration (h)') # + dict_unit[str(int(freq_factor))] + ')')
            if reg_todo is not 'CSA':
                subpl2.text(7, 2.0, reg_todo)
            else:
                subpl2.text(7, 2.0, 'SSZ')
            subpl2.set_xticks([1, 3, 6, 10, 16])
            subpl2.set_xticklabels([3, 9, 18, 30, 48])#['{:d}'.format(int(i)) for i in plotting_dur_bins[0:22:2]/3])
            # plt.colorbar(cms, shrink = 0.8)
            # subpl2.set_ylim(x_lim[0], x_lim[1])
            # subpl2.yaxis.get_major_ticks().label.set_fontsize(10)
            plt.ylim(1, 11)
            subpl2.set_yticks(range(1, 12, 1))
            ystr = []
            for i in plotting_intens_bins[1:12]:
                if i < 1:              
                    ystr.append('{:.1f}'                                                             .format(i))
                else:
                    ystr.append('{:d}'.format(int(i)))
            subpl2.set_yticklabels(ystr)

            if wet_or_dry_spells == 'wet':
                # subpl2.set_yscale('log')
                # subpl2.set_xlim(1, ylim)
                subpl2.set_xlim(0, ylim)
                # subpl2.set_xticks(range(3, ylim+3,3))
            plt.tick_params(axis='both', which='major', labelsize=10)
            if num_lines > 1:
                if len(numref_list) > 1:
                    if np.mod(numplot - 1, (numref_list[1] - 1)) == 0:
                        if integral_over_intensity:
                            subpl2.set_ylabel('number of dry spells/year')
                        else:
                            subpl2.set_ylabel('mean intensity (mm/' + dict_unit[str(int(freq_factor))] + ')')
                        if (numplot) > (num_lines - 1) * (numref_list[1] - 1):
                            subpl2.set_xlabel('duration (' + dict_unit[str(int(freq_factor))] + ')')
                else:
                    if numplot in [1, 3]:
                        subpl2.set_ylabel('cumulated rain (mm)')# + dict_unit[str(int(freq_factor))] + ')')
                    #if (numplot) > 2:
            elif (num_lines == 1) & (numplot == 1):
                if integral_over_intensity:
                    subpl2.set_ylabel('number of dry spells/year')
                else:
                    subpl2.set_ylabel('cumulated rain (mm)')
                    #subpl2.set_ylabel('mean intensity (mm/' + dict_unit[str(int(freq_factor))] + ')')

            if not integral_over_intensity:
                if len(numref_list) > 1:
                    cdt = (np.mod(numplot, numref_list[1] - 1) == 0) & (numplot == numref_list[1] - 1)
                else:
                    cdt = np.mod(numplot, 2) == 0
                if cdt:
                    if diff_btw_res:
                        colorbar_axes = plt.gcf().add_axes([0.915, 0.52, 0.02, 0.45])
                        colorbar = plt.colorbar(cms, colorbar_axes,)
                        #                        ticks=[-0.01, -0.001, -0.0001, -0.00001, 0.00001, 0.0001, 0.001,
                        #                               0.01], format='%.0e')
                        #colorbar.ax.tick_params(labelsize=8)
                    else:
                        colorbar_axes = plt.gcf().add_axes([0.915, 0.15, 0.02, 0.7])
                        colorbar = plt.colorbar(cms, colorbar_axes,)
                        colorbar.set_ticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10])
                        colorbar.set_ticklabels(['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5', '10'])#0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001])
                        #colorbar.ax.tick_params(labelsize=8)
                if len(numref_list) >= 1:
                    if numplot == numref_list[0] and (diff_btw_res is not None):
                        colorbar_axes = plt.gcf().add_axes([0.915, 0.02, 0.02, 0.45])
                        plt.colorbar(cms, colorbar_axes, )#ticks=[0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001])
                # colorbar_axes = plt.gcf().add_axes([0.93, 0.15, 0.02, 0.7])
                # else:
                #   colorbar_axes = plt.gcf().add_axes([0.95, 0.525, 0.45, 0.02])
                # colorbar = plt.colorbar(cms, colorbar_axes)
                # colorbar.ax.tick_params(labelsize=8)
                # colorbar.ax.xaxis.set_label_position('top')
            if integral_over_intensity:
                subpl2.legend(prop={'size': 16}, loc='best')
            print totnum
            if replace_names:
                if (diff_btw_res is not None) and (numplot not in numref_list):
                    subpl2.set_title('{}, {}\%'.format(resol[res], int(totnum_diff)))
                else:
                    subpl2.set_title('{}, {:.1f}'.format(resol[res], totnum_per_cell_per_year))
            else:
                if res in resol.keys():
                    res = resol[res]
                if integral_over_intensity:
                    subpl2.set_title('')  # {}'.format(season))  # , self.config['totnb'][res][0]))
                else:
                #     if (res_ref == 'TRMM'):
                #         if res == 'CP4A':
                #             res = 'b. CP4A-TRMM'
                #         elif res == 'N512r':
                #             res = 'c. R25-TRMM'
                #     else:
                #         if res == 'CP4A':
                #             res = 'e. CP4A-CMORPH'
                #         elif res == 'N512r':
                #             res = 'f. R25-CMORPH'
                    subpl2.set_title('{}, nb={:.1f}'.format(res,  # season.upper(),
                                                            totnum_per_cell_per_year),
                                      fontsize=14)  # , self.config['totnb'][res][0]))
        plt.tight_layout(rect=[0, 0, 0.9, 1], pad=0.05)
        plt.savefig(target, dpi=300)  # '2dhists_{}_{}_{}.png'.format(reg, season, '_'.join(r for r in res_todo)))


def generate_sig_mask(hist_ref, hist, res_ref, res, config, intensity_bins,
                                duration_bins, following_drought_bins=None):
    histo_dict = {res: hist, res_ref: hist_ref}
    cube_dict = {}
    simulated_dict = {}
    for reso in histo_dict:
        start_year = min([int(i) for i in histo_dict[reso].keys()])
        end_year = max([int(i) for i in histo_dict[reso].keys()])
        hist2plot = get_mean_10_year_change(histo_dict[reso].copy())
        cube_dict[reso] = hist2cube(hist2plot, intensity_bins, duration_bins, following_drought_bins)
        all_year_indices = generate_year_indices(start_year, end_year, config['bootstrap_samples_split'],
                          config['bootstrap_no_years'], config['bootstrap_filename'], config['master_ref_year'],
                          incomplete_years=False, bootstrap_year_set=[], print_diags=False)
        simulated_dict[reso] = generate_bootstrapped_cubes(histo_dict[reso].copy(), all_year_indices, intensity_bins,
                                                           duration_bins, following_drought_bins=None)

    sig_mask = generate_global_pvalue_mask(cube_dict[res_ref].copy(), simulated_dict[res_ref].copy(),
                                           cube_dict[res].copy(), simulated_dict[res].copy(), )

    return sig_mask


def test_event_change_significance(hist_ref, hist, res_ref, res,
                                   config, event, nb_points, nb_points_ref, season_length):
    histo_dict = {res: hist, res_ref: hist_ref}
    nb_points_dict = {res: nb_points, res_ref: nb_points_ref}
    cube_dict = {}
    simulated_dict = {}
    for reso in histo_dict:
        start_year = min([int(i) for i in histo_dict[reso].keys()])
        end_year = max([int(i) for i in histo_dict[reso].keys()])
        hist2plot = get_mean_10_year_change(histo_dict[reso].copy())
        _, event_nb = event.freq_or_nb_of_event(hist2plot,
                                                nb_points=nb_points_dict[reso],
                                                diff='percdiff', season_length=season_length)
        cube_dict[reso] = event2cube(event_nb)
        all_year_indices = generate_year_indices(start_year, end_year, config['bootstrap_samples_split'],
                          config['bootstrap_no_years'], config['bootstrap_filename'], config['master_ref_year'],
                          incomplete_years=False, bootstrap_year_set=[], print_diags=False)
        simulated_dict[reso] = generate_bootstrapped_event_cubes(histo_dict[reso].copy(), all_year_indices, event,
                                                                 nb_points_dict[reso], season_length, )
        # print 'simulated_dict[reso]', simulated_dict[reso]
    #print 'cube_dict[res_ref].data', cube_dict[res_ref].data
    p_val = estimate_local_pvalue_from_bootstrap(cube_dict[res_ref].copy(), simulated_dict[res_ref].copy(),
                                                 cube_dict[res].copy(), simulated_dict[res].copy(), )
    print 'p_val', p_val
    return p_val


### ME2: add this function
def get_histo_significance_interval(hist, res, config, bins, percentiles=[2.5, 97.5], bootstrap_or_centiles='bootstrap'):
    histo_dict = {res: hist, }
    cube_dict = {}
    simulated_dict = {}
    for reso in histo_dict:
        start_year = min([int(i) for i in histo_dict[reso].keys()])
        end_year = max([int(i) for i in histo_dict[reso].keys()])
        hist2plot = get_mean_10_year_change(histo_dict[reso].copy())
        cube_dict[reso] = hist2cube(hist2plot, bins)
        if bootstrap_or_centiles == 'bootstrap':
           all_year_indices = generate_year_indices(start_year, end_year, config['bootstrap_samples_split'],
                          config['bootstrap_no_years'], config['bootstrap_filename'], config['master_ref_year'],
                          incomplete_years=False, bootstrap_year_set=[], print_diags=False)
           simulated_dict[reso] = generate_bootstrapped_hist_cubes(histo_dict[reso].copy(), all_year_indices, bins)
        else:
           simulated_dict[reso] = generate_cube_from_hist_with_year_dim(histo_dict[reso].copy(), bins)
        # print 'simulated_dict[reso]', simulated_dict[reso]
    #print 'cube_dict[res].data', cube_dict[res].data
    if bootstrap_or_centiles == 'bootstrap':
       print 'centiles over bootstrap', res
       low_high_estimate = simulated_dict[res].collapsed('bootstrap_sample_number', iana.PERCENTILE, percent=percentiles)
    else:
       print 'centiles over the years/ens. spread', res
       print simulated_dict[res].coord('years')
       low_high_estimate = simulated_dict[res].collapsed('years', iana.PERCENTILE, percent=percentiles)
    mean = cube_dict[res]
    print 'low_high_estimate', low_high_estimate
    # p_val = estimate_local_pvalue_from_bootstrap(cube_dict[res_ref].copy(), simulated_dict[res_ref].copy(),
    #                                              cube_dict[res].copy(), simulated_dict[res].copy(), )
    # print 'p_val', p_val
    return low_high_estimate, mean

####update11/10
def test_hist_change_significance(hist_ref, hist, res_ref, res, config_ref, config, bins, ):
    histo_dict = {res: hist, res_ref: hist_ref}
    config_dict = {res: config, res_ref: config_ref}
    cube_dict = {}
    simulated_dict = {}
    for reso in histo_dict:
        start_year = min([int(i) for i in histo_dict[reso].keys()])
        end_year = max([int(i) for i in histo_dict[reso].keys()])
        hist2plot = get_mean_10_year_change(histo_dict[reso].copy())
        cube_dict[reso] = hist2cube(hist2plot, bins)
        config = config_dict[reso]
        all_year_indices = generate_year_indices(start_year, end_year, config['bootstrap_samples_split'],
                                                 config['bootstrap_no_years'], config['bootstrap_filename'],
                                                 config['master_ref_year'],
                                                 incomplete_years=False, bootstrap_year_set=[], print_diags=False)
        simulated_dict[reso] = generate_bootstrapped_hist_cubes(histo_dict[reso].copy(), all_year_indices, bins)
        # print 'simulated_dict[reso]', simulated_dict[reso]
    print 'cube_dict[res].data', cube_dict[res].data

    mean = cube_dict[res]
    #print 'low_high_estimate', low_high_estimate
    p_val = estimate_local_pvalue_from_bootstrap(cube_dict[res_ref].copy(), simulated_dict[res_ref].copy(),
                                                 cube_dict[res].copy(), simulated_dict[res].copy(), )
    p_val = np.ma.masked_where((cube_dict[res_ref].data <= 0.0005) & (cube_dict[res].data <= 0.0005), p_val)
    print 'p_val', p_val
    return p_val


####end update11/10
def get_mean_10_year_change(dict_by_years):
    min_year = int(min(dict_by_years.keys()))
    nb_year = len(dict_by_years.keys())
    mean_hist2plot = dict_by_years[str(min_year)].copy()
    for year in range(min_year+1, min_year+nb_year):
        mean_hist2plot += dict_by_years[str(year)]
    mean_hist2plot /= nb_year
    return mean_hist2plot


class EventsDefinitions:
    def __init__(self, name, duration_def=None, intensity_def=None, drought_def=None,
                 dim_dur=None, dim_int=None, dim_drought=None, shift=0.2, histoclass=None,
                 proportion_dry_days_in_season=False, units='nb events/year', hist_type='2d'):
        '''
        A class to define event types
        :param str name
        :param array/list duration_def: indices (start, end, step) of duration bins, end not included
        :param array/list intensity_def: same as above for intensity
        :param array/list drought_def: same as above for drought
        :param float shift: position of the event on the xaxis of figure
        :param int dim_dur: which dimension of the histogram is duration
        :param int dim_int: which dimension of the histogram is intensity
        :param int dim_drought: which dimension of the histogram is drought
        :param float shift: place on the figure (fraction of xaxis)
        :param PrecipExtreme histoclass:
        :param bool proportion_dry_days_in_season: for dry spells
        :param str units: default is 'nb events/year'. if 'nb in 10 years', then nb_years/10,
        :param str hist_type: if the histogram to fetch is from a '1d' calculation (only intensity bins) or '2d': intensity/duration/drought following event
        '''

        self.event_name = name
        if duration_def is not None:
            self.duration_def = slice(duration_def[0], duration_def[1], duration_def[2])
        self.intensity_def = slice(intensity_def[0], intensity_def[1], intensity_def[2])
        if drought_def is not None:
            self.drought_def = slice(drought_def[0], drought_def[1], drought_def[2])
        self.dim_duration = dim_dur
        self.dim_intensity = dim_int
        self.dim_drought = dim_drought
        self.shift = shift
        self.histo_class = histoclass
        self.proportion_dry_days_in_season = proportion_dry_days_in_season
        self.units = units
        self.hist_type = hist_type


    def freq_or_nb_of_event_per_hist(self, hist3d, subset, nb_points, season_length):
        if self.proportion_dry_days_in_season:
            nb_dry_days_per_bins = np.sum(hist3d[:, subset[1], subset[2]],
                                          axis=(1, 2)) * self.histo_class.duration_bins[1:]
            event_type_nb = np.sum(nb_dry_days_per_bins[subset[0]]) / (np.float(nb_points) * season_length)
        else:
            if len(subset) == 1:
                event_type_nb = np.sum(hist3d[subset[0],]) / np.float(nb_points)
            else:
                event_type_nb = np.sum(hist3d[subset[0], subset[1], subset[2]]) / np.float(nb_points)
        tot_events = np.sum(hist3d) / np.float(nb_points)
        event_type_freq = event_type_nb / tot_events

        return event_type_freq, event_type_nb


    def freq_or_nb_of_event(self, hist3d, hist3d_ref=None, nb_points=None, diff='diff', season_length=90):
        '''
        :param event_name:
        :param hist3d: ()
        :param hist3d: ()
        :param nb_points: (nb_years, nb_years_ref) normalisation factor, if nb_years_ref == None, assumes
        :return:
        '''
        hist_shape = np.array(hist3d).shape
        if len(hist_shape) == 3:
            subset = [0, 0, 0]
        elif len(hist_shape) == 1:
            subset = [0, ]
            if self.dim_intensity > 0:
                raise UserWarning('hist shape is 1, so intensity should be in dimension 0')
        else:
            raise Exception('hist shape {} not supported, only shape 1 or 3 supported'.format(len(hist_shape)))

        if self.dim_duration is not None:   subset[self.dim_duration] = self.duration_def
        if self.dim_intensity is not None:   subset[self.dim_intensity] = self.intensity_def
        if self.dim_drought is not None:   subset[self.dim_drought] = self.drought_def

        # print(subset)
        if nb_points[0]:
            event_type_freq, event_type_nb = self.freq_or_nb_of_event_per_hist(hist3d, subset,
                                                                               nb_points[0], season_length)
        else:
            raise UserWarning("please specify number of years or points to normalise")
        if hist3d_ref is not None:
            nb_year_ref = np.float(nb_points[0])
            if nb_points[1]:
                nb_year_ref = np.float(nb_points[1])
            event_type_freq_ref, event_type_nb_ref = self.freq_or_nb_of_event_per_hist(hist3d_ref, subset,
                                                                                  nb_year_ref, season_length)

            if 'percdiff' in diff:
                event_type_nb = (event_type_nb - event_type_nb_ref) / event_type_nb_ref * 100
                event_type_freq = (event_type_freq - event_type_freq_ref) / event_type_freq_ref * 100
            else:
                event_type_nb = event_type_nb - event_type_nb_ref
                event_type_freq = event_type_freq - event_type_freq_ref

        return event_type_freq, event_type_nb


def save_time_series_of_event_types(event_list, reg_todo, plot_type, colorlwdict, figdir, season=None, bigregion=None,
                     subregion_files=None, nb_freq_per_day=1., diff_btw_res=None, res_shift_dict=None, ylim_top=5.5,
                     xlim_top=1.6, no_map=False, bootstrap_config=None):
    '''

    :param event_list: a list of events defined by the class EventsDefinition
    :param reg_todo: list of region names to do, not used if subregion_files
    :param plot_type: 'diff' or 'percdiff if diff_btw_res
    :param colorlwdict: dictionary to define the color and line width/style of each resolution
    :param season: for title
    :param tuple bigregion: region limits, for set_extent
    :param subregion_files: if exist, then masked regions used from files
    :param int nb_freq_per_day:
    :param diff_btw_res: list of resolutions, ref = list[:-1]
    :param dict_names: if you want to change the standard resolution names
    :return:
    '''
    freqs = {'1': 'd', '8': '3h', '24': 'h', '96': '15mn'}
    r_letters = {'SEN': 'a. ', 'CSA': 'b. ', 'GOG': 'd. ', 'NSA': 'b.', 'GSL': 'c.'}
    print 'diff_btw_res', diff_btw_res, bool(diff_btw_res)

    if diff_btw_res:
        # plot the difference: diff_btw_res[1]-diff_btw_res[0], diff_btw_res[3]-diff_btw_res[2], ...
        res_todo = diff_btw_res[1::2]  # take the first and then every 2 datasets
        diff_res_todo = diff_btw_res[::2]  # take reference datasets every 2 datasets.
    else:
        res_todo = event_list[0].histo_class.resolutions
    hfiter = False

    for reg in reg_todo:
        if diff_btw_res:
            target = "{}/{}_event_time_series_{}_{}_{}_{}".format(figdir, reg, plot_type, season,
                                                       int(24. / nb_freq_per_day), diff_btw_res[0])
        else:
            target = "{}/{}_event_time_series_{}_{}_{}".format(figdir, reg, plot_type, season,
                                                    int(24. / nb_freq_per_day))

        #   if (len(reg_todo) == 1) & (len(frequencies)>1):
        ################ adddddd ############
        print target
        num_plots = len(event_list)
        from matplotlib import style
        style.use('fivethirtyeight')
        fig = plt.figure(figsize=(10, 3))
        start_plot = 1
        num_lines = 2
        numcol = 3
        # if reg_todo:
        #     for r in srex_reg.keys():
        #         x, y = srex_reg[r].exterior.xy
        #         subpl.plot(x, y, color='k')  # plot regional box lines
        #         xc, yc = np.array(srex_reg[r].exterior.centroid)
        #         subpl.text(xc, yc, r, horizontalalignment='center', verticalalignment='center', color='k',
        #                    fontdict={'weight': 'bold'})  # add region codes to boxes

        for regcount, event in enumerate(event_list):
            units = 'nb/season'
            numplot = start_plot + regcount
            subpl2 = fig.add_subplot(num_lines, numcol, numplot)

            for nb_r, r in enumerate(res_todo):
                if not ((event.event_name == 'Intense\ndays\n(OND)') and (reg == 'NSA')):
                    if not res_shift_dict:
                        res_shift = np.float(nb_r) * 0.12 / (len(res_todo))
                    else:
                        res_shift = res_shift_dict[r]
                    print event.event_name
                    hist2plot_all, nb_years, nb_spatial_pts, _ = event.histo_class.get_hist2plot(r, reg, hfiter,
                                                                                                 d1_or_2d=event.hist_type)
                    if event.hist_type == '2d':
                        hist2plot_all_years = hist2plot_all[1].copy()
                    else:
                        hist2plot_all_years = hist2plot_all[0].copy()
                    year_list = []
                    event_list_by_year = []
                    for year in sorted(hist2plot_all_years):
                        print year
                        year_list.append(year)
                        barwidth = 0.12 / (len(res_todo) + 2)  # -2
                        season_length = 90.
                        if r in ['TRMM', 'CMORPH']:
                            barwidth = 0.01
                            season_length = 92.
                        # if 'floods' in event.event_name:
                        scaling_factor = 1.
                        if event.units == 'nb in 10 years':
                            scaling_factor = 0.1
                            units = 'nb in 10 seasons'
                        # print nb_years, nb_spatial_pts
                        event_type_freq, event_type_nb = event.freq_or_nb_of_event(hist2plot_all_years[year],
                                                                                   nb_points=(nb_spatial_pts*scaling_factor, None),
                                                                                   season_length=season_length)
                        event_list_by_year.append(event_type_nb)
                    print year_list
                    print event_list_by_year
                    subpl2.plot(year_list, event_list_by_year, color=colorlwdict[r][0], alpha=1.0,
                               label=r,)
                    title = ' '.join(event.event_name.split('\n'))
                    if 'Dry\nspells' in event.event_name:
                        units = 'fraction of dry spell days in season'
                    plt.title('{} ({})'.format(title, units), fontsize=8)
                    plt.xticks(rotation=70)
                    if regcount < 3:
                        subpl2.set_xticklabels([])
        plt.savefig(target + '.pdf', bbox_inches='tight')  # dpi=300, bbox_inches='tight')
        plt.savefig(target + '.png', dpi=300, bbox_inches='tight')


def save_event_types(event_list, reg_todo, plot_type, colorlwdict, figdir, season=None, bigregion=None,
                     subregion_files=None, nb_freq_per_day=1., diff_btw_res=None, res_shift_dict=None, ylim_top=5.5,
                     xlim_top=1.6, no_map=False, bootstrap_config=None):
    '''

    :param event_list: a list of events defined by the class EventsDefinition
    :param reg_todo: list of region names to do, not used if subregion_files
    :param plot_type: 'diff' or 'percdiff if diff_btw_res
    :param colorlwdict: dictionary to define the color and line width/style of each resolution
    :param season: for title
    :param tuple bigregion: region limits, for set_extent
    :param subregion_files: if exist, then masked regions used from files
    :param int nb_freq_per_day:
    :param diff_btw_res: list of resolutions, ref = list[:-1]
    :param dict_names: if you want to change the standard resolution names
    :return:
    '''
    freqs = {'1': 'd', '8': '3h', '24': 'h', '96': '15mn'}
    r_letters = {'SEN': 'a. ', 'CSA': 'b. ', 'GOG': 'd. ', 'NSA': 'b.', 'GSL': 'c.'}
    print 'diff_btw_res', diff_btw_res, bool(diff_btw_res)
    if diff_btw_res:
        # plot the difference: diff_btw_res[1]-diff_btw_res[0], diff_btw_res[3]-diff_btw_res[2], ...
        res_todo = diff_btw_res[1::2]  # take the first and then every 2 datasets
        diff_res_todo = diff_btw_res[::2]  # take reference datasets every 2 datasets.
    else:
        res_todo = event_list[0].histo_class.resolutions
    hfiter = False
    if diff_btw_res:
        target = "{}/{}_events_{}_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season,
                                                   int(24. / nb_freq_per_day), diff_btw_res[0])
    else:
        target = "{}/{}_events_{}_{}_{}".format(figdir, '_'.join(s for s in reg_todo), plot_type, season,
                                                int(24. / nb_freq_per_day))
    regions = reg_todo
    if reg_todo == 'masked':
        hfiter = True
        reg_todo = event_list[0].histo_class.hists.keys[0]['masked'].keys()
        regions = None
        target = "{}/{}_{}_{}_{}_{}_{}".format(figdir, '_'.join(s[4:-5] for s in reg_todo), res_todo[0], plot_type,
                                               nb_freq_per_day, freqs[str(int(nb_freq_per_day))], season)
    #   if (len(reg_todo) == 1) & (len(frequencies)>1):
    ################ adddddd ############
    print target
    num_plots = len(reg_todo)
    from matplotlib import style
    style.use('fivethirtyeight')
    fig, subpl, numcol, start_plot, num_lines, col_order = set_up_fig(num_plots, reg_todo=regions, bigregion=bigregion,
                                                                      subregion_files=subregion_files, no_map=no_map)
    # if reg_todo:
    #     for r in srex_reg.keys():
    #         x, y = srex_reg[r].exterior.xy
    #         subpl.plot(x, y, color='k')  # plot regional box lines
    #         xc, yc = np.array(srex_reg[r].exterior.centroid)
    #         subpl.text(xc, yc, r, horizontalalignment='center', verticalalignment='center', color='k',
    #                    fontdict={'weight': 'bold'})  # add region codes to boxes

    if no_map == False:
        if regions:
            if len(regions) > 1:
                subpl.set_extent(bigregion, ccrs.PlateCarree())  # Africa
    for regcount, reg in enumerate(reg_todo):
        numplot = start_plot + regcount
        subpl2 = fig.add_subplot(num_lines, numcol, numplot)
        if season:
            season = reg
        print reg
        xticks = []
        xticklabels = []
        if reg is not None:
            plt.xlim((0.07, xlim_top))
            if reg == 'GSL':
                plt.ylim((0, 28))
            elif reg == 'NSA':
                plt.ylim((0, 10))
            else:
                plt.ylim((0, ylim_top))
            subpl2.spines['right'].set_visible(False)
            subpl2.spines['top'].set_visible(False)

            bin_oc_tot = []
            arrow_dict = {}
            for nb_r, r in enumerate(res_todo):
                for event in event_list:
                    if not ((event.event_name == 'Intense\ndays\n(OND)') and (reg == 'NSA')):
                        if not res_shift_dict:
                            res_shift = np.float(nb_r) * 0.12 / (len(res_todo))
                        else:
                            res_shift = res_shift_dict[r]
                        print event.event_name
                        if nb_r == 0:
                            xticks = np.append(xticks, event.shift + 0.)
                            xticklabels = np.append(xticklabels, event.event_name)
                        hist2plot_all, nb_years, nb_spatial_pts, _ = event.histo_class.get_hist2plot(r, reg, hfiter, d1_or_2d=event.hist_type)
                        if event.hist_type == '2d':
                            hist2plot_all_years = hist2plot_all[1].copy()
                        else:
                            hist2plot_all_years = hist2plot_all[0].copy()

                        hist2plot = get_mean_10_year_change(hist2plot_all_years.copy())
                        barwidth = 0.12 / (len(res_todo)+2) #-2
                        season_length = 90.
                        if r in ['TRMM', 'CMORPH']:
                            barwidth = 0.01
                            season_length = 92.
                        # if 'floods' in event.event_name:
                        if event.units == 'nb in 10 years':
                            nb_years = nb_years / 10.
                        print nb_years, nb_spatial_pts
                        event_type_freq, event_type_nb = event.freq_or_nb_of_event(hist2plot,
                                                                                   nb_points=(nb_years * nb_spatial_pts, None),
                                                                                   diff='percdiff',
                                                                                   season_length=season_length)
                        # subpl2.bar(event.shift + res_shift, event., color=colorlwdict[r][0], alpha=0.8)
                        subpl2.bar(event.shift + res_shift, event_type_nb, color=colorlwdict[r][0], alpha=1.0,
                                   label=r, width=barwidth, zorder=10)

                        if (diff_btw_res is not None) & (r not in ['TRMM', 'CMORPH']):
                            res_shift = res_shift_dict[diff_res_todo[nb_r]]
                            hist2plot_ref_all, nb_years_ref, nb_spatial_pts_ref, _ = event.histo_class.get_hist2plot(
                                diff_res_todo[nb_r],
                                reg, hfiter, d1_or_2d=event.hist_type)
                            # if 'floods' in event.event_name:
                            if event.units == 'nb in 10 years':
                                nb_years_ref = nb_years_ref / 10.
                            if event.hist_type == '2d':
                                hist2plot_ref_all_years = hist2plot_ref_all[1]
                            else:
                                hist2plot_ref_all_years = hist2plot_ref_all[0]
                            hist2plot_ref = get_mean_10_year_change(hist2plot_ref_all_years.copy())
                            print('calculating significance')

                            fontweight = 'normal'
                            if bootstrap_config is not None:
                                significance_of_change = test_event_change_significance(hist2plot_ref_all_years,
                                                                                    hist2plot_all_years,
                                                                                    res='future',
                                                                                    res_ref='present',
                                                                                    config=bootstrap_config,
                                                                                    event=event,
                                                                                    nb_points=(nb_years * nb_spatial_pts, None),
                                                                                    nb_points_ref=(nb_years_ref * nb_spatial_pts_ref, None),
                                                                                    season_length=season_length)
                                if significance_of_change[0] < 0.01:
                                    fontweight = 'bold'

                            print nb_years_ref, nb_spatial_pts_ref
                            #print 'histoplot_ref', hist2plot_ref
                            _, nb_of_event_ref = event.freq_or_nb_of_event(hist2plot_ref, nb_points=(nb_years_ref * nb_spatial_pts_ref, None),
                                                      season_length=season_length)

                            _, event_type_change = event.freq_or_nb_of_event(hist2plot, hist2plot_ref,
                                                                             (nb_years * nb_spatial_pts,
                                                                              nb_years_ref * nb_spatial_pts_ref),
                                                                             diff=plot_type,
                                                                             season_length=season_length)
                            format1 = '{}{:.1f}'
                            if 'percdiff' in plot_type:
                                format1 = '{}{:.0f}%'
                            sign = ''
                            zorder = 9
                            if event_type_change > 0:
                                sign = '+'
                                zorder = 11
                            subpl2.bar(event.shift + res_shift, nb_of_event_ref,
                                           color=colorlwdict[diff_res_todo[nb_r]][0], alpha=1.0,
                                           label=r, width=barwidth, zorder=zorder)
                            text_shift = 0.07
                            if ylim_top < 3:
                                text_shift = 0.07/10.

                            subpl2.text(event.shift + res_shift-0.02, event_type_nb + text_shift,
                                        format1.format(sign, event_type_change), color='k', fontsize=8, zorder=12,
                                        weight=fontweight)
            subpl2.set_xticks(xticks)
            subpl2.set_xticklabels(xticklabels, fontsize=7)  # , rotation=15)
        # ylabel = 'event/year/grid point'
        # if (np.mod(numplot - 1, numcol) == 0):
        #     plt.tick_params(
        #         axis='y',  # changes apply to the y-axis
        #         labelleft='on')
        #  #   subpl2.set_ylabel(ylabel)
        # elif (np.mod(numplot, numcol) == 0):
        #     plt.tick_params(
        #         axis='y',  # changes apply to the y-axis
        #         labelleft='off',
        #         labelright='on')
        #  #   subpl2.set_ylabel(ylabel)
        # else:
        #     plt.tick_params(
        #         axis='y',  # changes apply to the y-axis
        #         labelleft='on')  # 'off')
        #     #subpl2.set_ylabel(ylabel)

        if reg in r_letters:
            if reg == 'CSA':
                reg = r_letters[reg] + 'SSZ'
            else:
                reg = r_letters[reg] + reg
        else:
            reg = ''
        # if 'GSL' in reg:
        #     plt.title(reg, y=0.95)
        # else:
        plt.title(reg, x=0.55, y=0.83, loc='center')  # 0.89 # # 'b. using CMORPH values' #'a. using TRMM values' #'b. using CMORPH values'
        # if regcount == 0:
        #     subpl2.legend(prop={'size': 6}, loc='best')
    plt.show()
    plt.savefig(target + '.pdf', bbox_inches='tight')  # dpi=300, bbox_inches='tight')
    plt.savefig(target + '.png', dpi=300, bbox_inches='tight')

def check_precip_mm_or_kg_m2_s1(cube):
    '''
    What it does  : checks units are mm, mm day-1 or kg.m-2.s-1
    Modified      :
    :param str freq: day or 3h for example  
    '''

    if (cube.units in ['kg m-2 s-1', 'mm day-1', 'mm h-1', 'mm (h)-1', 'mm (3h)-1']):
        outcube = cube
    elif (cube.units in ['mm']):
        if cube.attributes.has_key('short_name'):
            if cube.attributes['short_name'] == 'daily_rainfall':
                outcube = cube
                outcube.units = 'mm day-1'
        elif cube.attributes.has_key('comment'):
            if 'Disaggreg' in cube.attributes['comment']:
                outcube = cube
                outcube.units = 'mm h-1'
        else:
            outcube = cube
    elif (cube.units in ['10^-1 mm']):
        if any(m.method == 'sum' for m in cube.cell_methods):
            cube = cube * 10
            cube.units = 'mm day-1'
        else:
            cube = cube * 10
            cube.units = 'mm h-1'
    elif cube.units in ['kg m-2']:
        if cube.attributes.has_key('source'):
            if ('Met Office' in cube.attributes['source']) | \
                    ('cosmo_4_19_GPU_prototype' in cube.attributes['source']):
                outcube = cube
                outcube.units = 'mm h-1'  # nimrod
        else:
            rho = iris.coords.AuxCoord(1000.0, long_name='ref density', units='kg m-3')
            outcube = cube / rho
            outcube.convert_units('mm')
            outcube.units = 'mm day-1'

    elif cube.units in ['mm/3h', 'mm/3hr']:
        outcube = cube
        outcube.units = 'mm (3h)-1'
    else:
        raise ValueError('cube not in mm day-1, mm, mm h-1, kg/m2/s or kg m-2 but in {}'.format(cube.units))

    print 'outcube after units checked', outcube
    return outcube


def mask_cube_where(cube, cdt):
    '''Returns a cube with masked data where the condition cdt is met'''
    masked_cube = cube.copy()
    if type(masked_cube.data) == np.ndarray:
        masked_cube.data = np.ma.MaskedArray(cube.data, mask=False)
    masked_cube.data = np.ma.masked_where(cdt, masked_cube.data, copy=True)
    return masked_cube


def loadpp_regridded2eobs(infiledir, time_con):
    '''
   Loads rain and snow from the convection-permitting simulations
   and adds them up in total precipitation

   :param list infiledir: a list of length 1 with the directory where files are
   :param iris.Constraint time_con: time constraint for Giorgia's CPM simulations
    '''

    if 'gfosser' in infiledir[0]:
        precip = iris.load_cube(infiledir[0], time_con)
        return precip
    elif (len(glob.glob(infiledir[0] + '*_*_reg_rain.*')) > 0):
        rain = iris.load(infiledir[0] + '*_reg_rain.*', callback=callback_overwrite)
        snow = iris.load(infiledir[0] + '*_reg_snow.*', callback=callback_overwrite)
    else:
        rain = iris.load(infiledir[0] + '*.nc', 'stratiform_rainfall_rate', callback=callback_overwrite)
        snow = iris.load(infiledir[0] + '*.nc', 'stratiform_snowfall_rate', callback=callback_overwrite)

    rain = rain.concatenate_cube()
    snow = snow.concatenate_cube()
    for cube in [rain, snow]:
        if 'STASH' in cube.attributes.keys():
            cube.attributes['STASH'] = 'overwritten'
    cubelist = iris.cube.CubeList([rain, snow])
    iris.util.unify_time_units(cubelist)
    rain = cubelist[0]
    snow = cubelist[1]
    print rain
    print snow
    precip = rain + snow
    precip.long_name = 'precipitation'
    if ('alps' in infiledir[0]) & ('n512' not in infiledir[0].lower()):
        add_laea.load_laea(precip, '/project/hires_rcm/hadek/ALPS/RapdD_al05.etrs.laea_1989.nc')
    return precip


def set_up_fig(num_plots, reg_todo=None, bigregion=None, subregion_files=None, no_map=False, numref_list=None, season=None):
    col_order = None
    if (num_plots < 2):
        fig = plt.figure(figsize=(6, 7.2))  # (12,6))
        plt.subplots_adjust(wspace=0.2, hspace=0.2)
        num_lines = 1
        numcol = 1
        totcol = 1
        start_plot = 1
        subpl = None
    elif no_map:
        if (numref_list is not None):
            fig = plt.figure(figsize=(10, 2.5))
            if len(numref_list) > 0:
                num_lines = len(numref_list)
                totcol = 4  # 2#3
                num_lines = 1
            else:
                num_lines = 2
                totcol = 2
        else:
            fig = plt.figure(figsize=(10, 3))
            num_lines = 1 #num_plots//2+1
            totcol = 2 #num_plots//num_lines
        start_plot = 1
        subpl = None
    else:
        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(wspace=0.2, hspace=0.2)  # wspace=0.3,hspace=0.35)
        if num_plots <= 3:
            totcol = 2
        else:
            totcol = 3
        if np.mod(num_plots, 2) == 0:
            numcol = totcol - 1
            start_plot = 3
        else:
            numcol = totcol
            start_plot = 2
        if num_plots < 6:
            num_lines = 2
        else:
            raise Exception('can''t handle more than 6 plots')

        if no_map == False:
            subpl = fig.add_subplot(num_lines, numcol, 1, projection=ccrs.PlateCarree(), frameon=False)
            if bigregion == (-20, 12, 0, 25):
                mean_cube = iris.load_cube(
                    '/project/hires_rcm/sberthou/TRMM/PROCESSED/TRMM_Sahel_d_precip_JAS_mean_1_1998_2008.nc')
                cmap = cm.get_cmap('terrain_r', 9)
                levels = [1] + [i for i in range(2, 16, 2)]
                trmm_mask = iris.load_cube('/project/hires_rcm/CP4/TRMM/trmm.landmask.AFRICA.nc')
                trmm_mask = trmm_mask.regrid(mean_cube, iana.Nearest())
                mean_cube = mask_cube_where(mean_cube, trmm_mask.data < 0.5)
                cms = iplt.contourf(mean_cube, levels=levels, cmap=cmap, extend='both')
                plt.title('a. mean (TRMM)', weight='bold')
                colorbar_axes = fig.add_axes([0.07, 0.54, 0.013, 0.35])
                cbar = plt.colorbar(cms, colorbar_axes)
                cbar.set_label('mm day$^{-1}$', labelpad=-5, y=1.2, rotation=0)
                plt.tick_params(
                    axis='y',  # changes apply to the y-axis
                    labelleft='on',
                    labelright='off')  # 'off')
                subpl.outline_patch.set_visible(False)
            elif bigregion=='future_ESM':
                present_mean = load_ESM_mean('ESM-', season.upper())
                future_mean = load_ESM_mean('ESM-r', season.upper())
                cube2plot = future_mean - present_mean #)/present_mean*100
                cmap1 = cm.get_cmap('BrBG')
                cms = iplt.pcolormesh(cube2plot, vmin=-4, vmax=4, cmap=cmap1)
                iplt.contour(present_mean, levels=range(5, 20, 10), colors='k')
                colorbar_axes = fig.add_axes([0.07, 0.54, 0.013, 0.35])
                cbar = plt.colorbar(cms, colorbar_axes)
                plt.title('mm/day')

            # subpl.stock_img()
            subpl.coastlines('50m')
            # if not reg_todo:
            subpl.add_feature(cfeature.BORDERS, linewidth=0.5)
            # gl = subpl.gridlines(draw_labels=True)
            # gl.xlabels_top = False
            # gl.ylabels_right = False
            for ax_pos in ['top', 'bottom', 'left', 'right']:
                subpl.spines[ax_pos].set_visible(False)
            if not reg_todo:
                # pass
                ####update11/10
                col_order = plot_masked_subregions(subregion_files, subpl, season)
            else:
                for r in reg_todo:
                    x, y = srex_reg[r].exterior.xy
                    subpl.plot(x, y, color='k')  # plot regional box lines
                    xc, yc = np.array(srex_reg[r].exterior.centroid, )
                    if r == 'CSA':
                        r = 'SSZ'
                    subpl.text(xc, yc, r, horizontalalignment='center',
                               verticalalignment='center', color='red',
                           fontdict={'weight': 'bold'})  # add region codes to boxes
            if bigregion:
                if bigregion=='future_ESM':
                    polygon_list = []
                    for r in reg_todo:
                        polygon_list.append(srex_reg[r])
                    lonmin = min([min(r.exterior.xy[0]) for r in polygon_list])
                    lonmax = max([max(r.exterior.xy[0]) for r in polygon_list])
                    latmin = min([min(r.exterior.xy[1]) for r in polygon_list])
                    latmax = max([max(r.exterior.xy[1]) for r in polygon_list])
                    subpl.set_extent((lonmin-2, lonmax+2, latmin-2, latmax+2), ccrs.PlateCarree())
                else:
                    subpl.set_extent(bigregion, ccrs.PlateCarree())
    return fig, subpl, totcol, start_plot, num_lines, col_order

def load_ESM_mean(esm_name, season):
    '''
    loads mean rainfall for the 5 members
    :param esm_name: can be ESM- or ESM-r
    :return: mean of the 5 members
    '''

    files2load = '/project/hires_rcm/sberthou/{}?/PROCESSED/{}?_GLOBAL_d_precip_{}_mean_1_*_None.nc'.format(esm_name,
                                                                                                                    esm_name,
                                                                                                                    season)
    means = iris.load(files2load, callback=callback_overwrite)
    for m in range(5):
        means[m].add_aux_coord(iris.coords.AuxCoord(m + 1, long_name='member'))
        means[m] = iris.util.new_axis(means[m], 'member')
    all_means = means.concatenate_cube()
    mean = all_means.collapsed('member', iana.MEAN)
    if mean.units == 'kg m-2 s-1':
        mean = mean * 86400
        mean.units = 'mm day-1'

    return mean

####update11/10
def plot_masked_subregions(subregion_files, subpl, season=None):
    '''
    Plots the country dataset subdivided into regions according the the mask files in subregion_files
    '''
    cty = iris.load_cube(subregion_files.values()[0])
    if not isinstance(cty.data, np.ma.MaskedArray):
        masked_cube = cty.copy()
        masked_cube.data = np.ma.MaskedArray(cty.data, mask=False)
        cty = masked_cube
        # if ('alps' in subregion_files.values()[0]) & ('n512' not in subregion_files.values()[0]):
    #   add_laea.load_laea(cty, '/project/hires_rcm/hadek/ALPS/RapdD_al05.etrs.laea_1989.nc')
    tot_nb_reg = len(subregion_files.keys())
#    if 'IP' in subregion_files:
#       colors = plt.get_cmap('brewer_Set1_09')(range(0, tot_nb_reg))#'brewer_Set1_09')(range(0, tot_nb_reg))'Greys'
#    else:
    colors = plt.get_cmap('Greys', 10)(range(2, tot_nb_reg+2))#'brewer_Set1_09')(range(0, tot_nb_reg))'Greys'
    contour_levels = range(tot_nb_reg + 1)
    cmap, norm = from_levels_and_colors(contour_levels, colors)
    color_iter = iter(cmap(np.linspace(0, 1, tot_nb_reg)))
    col_order = {}
    for nbreg, subreg in enumerate(subregion_files.keys()):
        mcube = iris.load_cube(subregion_files[subreg])
        if not isinstance(mcube.data, np.ma.MaskedArray):
            masked_cube = mcube.copy()
            masked_cube.data = np.ma.MaskedArray(mcube.data, mask=False)
            mcube = masked_cube
        # if ('alps' in subregion_files[subreg].lower()) & ('n512' not in subregion_files[subreg]):
        #    add_laea.load_laea(mcube, '/project/hires_rcm/hadek/ALPS/RapdD_al05.etrs.laea_1989.nc')
        print subregion_files.values()[0], subreg
        mcube.data[mcube.data.mask==False] = nbreg
        col_order[subreg] = next(color_iter)
        iplt.pcolormesh(mcube, norm=norm, cmap=cmap)   

    #     colors = np.concatenate([plt.get_cmap('brewer_RdBu_11'  )(range(1, 10)),
    #              plt.get_cmap('brewer_YlOrBr_09')(range(5, 0, -1))])
    print 'plot_masked_subregions'
    print 'col_order', col_order
    # cmap = cm.get_cmap('Paired_r')
    #iplt.pcolormesh(cty, norm=norm, cmap=cmap)
    if season:
        plt.title(season.upper())
    for ax_pos in ['top', 'bottom', 'left', 'right']:
        subpl.spines[ax_pos].set_visible(False)

    return col_order
####end update11/10


def extract_region(cube, reg_lon_ce, reg_lat_ce):
    '''
    extracts the region given by reg_lon_ce coord_extend()
    '''
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    r_ppn = cube.intersection(reg_lat_ce).intersection(reg_lon_ce)
    return r_ppn


def hourly2daily(hourly_cube):
    '''
    Aggregates the hourly precipitation cube into daily_cube
    :param iris.cube: hourly cube to be aggregated in daily
    '''
    icat.add_day_of_year(hourly_cube, 'time')
    if (repr(hourly_cube.units) in ["Unit('kg/m2/s')", "Unit('kg m-2 s-1')"]):
        daily_data = hourly_cube.aggregated_by('day_of_year', iana.MEAN)
    elif (repr(hourly_cube.units) in ["Unit('mm')", "Unit('mm/h')"]):
        print 'summing hourly'
        daily_data = hourly_cube.aggregated_by('day_of_year', iana.SUM)
        daily_data.units = 'mm/day'
    else:
        raise (UserWarning)
    return daily_data


def plot_regions():
    srex_reg = {'WSA': Polygon(((-7., 6), (-17.0, 6.00), (-17., 13.), (-7., 13.))),
                'CSA': Polygon(((6., 9.0), (-7.0, 9.00), (-7., 13.), (6., 13.))),
                'GOG': Polygon(((6., 4.0), (-7.0, 4.00), (-7., 9.), (6., 9.))),
                'NSA': Polygon(((10., 13.0), (-17.0, 13.0), (-17., 19.), (10., 19.))),
                'BIG': Polygon(((40., 8.0), (-17.0, 8.00), (-17., 19.), (40., 19.))),
                'B': Polygon(((2.975, 8.525), (1.025, 8.525), (1.025, 10.475), (2.975, 10.475))),
                'M': Polygon(((-1.975, 14.775), (-1.025, 14.775), (-1.025, 15.725), (-1.975, 15.725))),
                'N': Polygon(((3.225, 12.775), (1.775, 12.775), (1.775, 14.225), (3.225, 14.225))),
                'JOS': Polygon(((10., 9.0), (6.0, 9.0), (6.0, 13.), (10., 13.))), }

    fig = plt.figure()
    subpl = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    subpl.stock_img()
    subpl.coastlines('50m')
    iplt.contour(orog, levels=[500, 1000, 2000, 3000, 4000], colors='k', linewidths=0.5)
    for r in srex_reg.keys():
        x, y = srex_reg[r].exterior.xy
        subpl.plot(x, y, color='k')
        xc, yc = np.array(srex_reg[r].exterior.centroid)
        subpl.text(xc, yc, r, horizontalalignment='center', verticalalignment='center', color='k',
                   fontdict={'weight': 'bold'})
    subpl.set_extent((-20, 50, 0, 30), ccrs.PlateCarree())
    plt.savefig('/data/users/sberthou/PRECIP_N512/Sahel/WAF_regions.png', dpi=300)
    plt.show()


def plot_percent_change_between_centiles(ax, hist_to_plot_fin_ref, hist_to_plot_fin, bins, between, ylim_max, color='k',
                                         shift=0, plot_only='all', ylim_top=0.025):
    '''
    
    
    :param ax: 
    :param hist_to_plot_fin_ref: 
    :param hist_to_plot_fin: 
    :param bins: 
    :param between: 
    :param color: 
    :return: 
    '''

    percent_change = hist_to_plot_fin / hist_to_plot_fin_ref * 100
    lower_change = np.sum(hist_to_plot_fin[bins < between[0]]) / np.sum(hist_to_plot_fin_ref[bins < between[0]]) * 100
    middle_change = np.sum(hist_to_plot_fin[(between[0] <= bins) & (bins < between[1])]) / np.sum(
        hist_to_plot_fin_ref[(between[0] <= bins) & (bins < between[1])]) * 100
    upper_change = np.sum(hist_to_plot_fin[bins >= between[1]]) / np.sum(hist_to_plot_fin_ref[bins >= between[1]]) * 100
    if plot_only == 'upper':
        future_hist = hist_to_plot_fin_ref + hist_to_plot_fin
        print('current contrib of upper 99.99:',
              np.sum(hist_to_plot_fin_ref[bins >= between[1]])/np.sum(hist_to_plot_fin_ref)*100)
        print('future contrib of upper 99.99:',
              np.sum(future_hist[bins >= between[1]])/np.sum(future_hist)*100)
    # subpl.arrow(between[0], -0.0028+ shift * 0.5, 0, , head_width=0.05, head_length=0.1, fc='k', ec='k')
    # plt.axvline(between[0], color=color, linewidth=2, linestyle='-.')
    # plt.text(between[0] * 1.05, -0.0028 + shift * 0.5, 'p99', color=color, fontsize=6)
    # plt.axvline(between[1], color=color, linewidth=2, linestyle='--')
    # plt.text(between[1] * 1.05, -0.0028 + shift * 0.5, 'p99.99', color=color, fontsize=6)

    # sign = value / np.abs(value)
    #             add_sign = {-1: '', 1: '+'}
    #             ax.text(annotx+np.random.rand()*2, sign*(annoty+np.random.rand()*0.0002),
    #                     '{}{:d}\%'.format(add_sign[sign], int(value)), fontsize=10,
    #
    #                    color=color,)
    ypos = -ylim_top/0.23
    # ypos = -0.0026 # -0.00195
    ypos_text =  ylim_top*0.01 # ypos - ypos/33
    # ypos_text = ypos - 0.00018
    # if not plot_only == 'upper':
    #     plt.plot([0, between[0] * 0.82], [ypos + shift, ypos + shift], lw=10, color=color, alpha=0.5)
    #     plt.plot([between[0] * 1.2, between[1] * 0.82], [ypos + shift, ypos + shift], lw=10, color=color, alpha=0.5)
    # plt.plot([between[1] * 1.2, ylim_max], [ypos + shift, ypos + shift], lw=10, color=color, alpha=0.5)

    if not plot_only == 'upper':
        if np.isfinite(lower_change):
            pos = between[0] / 2. * 0.4
            print(pos)
            plt.text(pos, ypos_text + shift, '{}{:d}%'.format('+' if lower_change > 0 else '', int(lower_change)),
                     fontsize=8)  # , color=color)
        if np.isfinite(middle_change):
            pos = (between[0] + between[1]) / 2. * 0.7
            print(pos)
            plt.text(pos, ypos_text + shift, '{}{:d}%'.format('+' if middle_change > 0 else '', int(middle_change)),
                     fontsize=8)  # , color=color)
    if np.isfinite(upper_change):
        #if plot_only == 'upper':
        #    pos = 30
        #else:
        pos = (between[1] + ylim_max) / 2. * 0.2
        print(pos)
        plt.text(pos, ypos_text + shift, '{}{:d}%'.format('+' if upper_change > 0 else '', int(upper_change)),
                 fontsize=8)  # , color=color)


def plot_one_third_two_thirds_vlines(ax, hist_to_plot_fin, bins, color='k',):

    hist_to_plot_fin_cum = np.cumsum(hist_to_plot_fin)/ hist_to_plot_fin.sum()
    below_one_third = np.where(hist_to_plot_fin_cum < 0.333)
    below_two_thirds = np.where(hist_to_plot_fin_cum > 0.666)
    one_third = below_one_third[0][-1]
    two_thirds = below_two_thirds[0][0]
    print one_third, bins[one_third]
    print two_thirds, bins[two_thirds]
    ax.arrow(bins[one_third], 0, 0, hist_to_plot_fin[one_third], head_width=0, head_length=0, fc=color, ec=color)
    ax.arrow(bins[two_thirds], 0, 0, hist_to_plot_fin[two_thirds], head_width=0, head_length=0, fc=color, ec=color)

    return bins[one_third], bins[two_thirds]

def plot_fraction_and_break_pts(ax, hist_to_plot_fin_ref, hist_to_plot_fin, bins, color='k', between=None):
    '''
     Adds on the x axis the values when differences reach 0 and add the fraction of over/underestimation integrated between these points
    :param ax : figure axis (ax = plt.add_subplot())
    '''
    percent_change = hist_to_plot_fin / hist_to_plot_fin_ref * 100
    where_zero = hist_to_plot_fin[:-1] * hist_to_plot_fin[1:]
    where_zero_are_init = np.where(where_zero <= 0)
    if where_zero_are_init[0][0] < 10:
        where_zero_are_init[0][0] = 0
    where_zero_are = [where_zero_are_init[0][0], ]
    print where_zero_are_init

    def consecutive(data, stepsize=1):
        return np.split(data, np.where((np.diff(data) > stepsize))[0] + 1)

    all_zeros = consecutive(where_zero_are_init[0], 2)
    print 'all_zeros', all_zeros
    where_zero_are = [a[0] for a in all_zeros]  # only take a place if not consecutive
    if all_zeros[-1][-1] < len(bins) - 3:
        where_zero_are.append(len(bins) - 3)

    if len(where_zero_are_init[0]) > 0:
        def consecutive(data, stepsize=1):
            return np.split(data, np.where((np.diff(data) > stepsize))[0] + 1)

        all_zeros = consecutive(where_zero_are_init[0], 10)
        print 'all_zeros', all_zeros
        where_zero_are = [a[0] for a in all_zeros]  # only take a place if not consecutive
        if all_zeros[-1][-1] < len(bins) - 3:
            where_zero_are.append(len(bins) - 3)

        print where_zero_are
        where_zero_are = (where_zero_are,)
        where_add_ticks = bins[where_zero_are]
        print where_add_ticks
        #        ax.set_xticks(where_add_ticks, minor=True)
        #         ax.set_xticklabels(['{:.1f}'.format(a) for a in where_add_ticks], minor=True, color='r', va='top')
        # integrated_percent_change = [percent_change[:where_zero_are[0]].sum(),]
        if len(where_zero_are[0]) > 1:
            integrated_percent_change = [
                hist_to_plot_fin[am1:am + 1].sum() / hist_to_plot_fin_ref[am1:am + 1].sum() * 100 for am1, am in
                zip(where_zero_are[0][:-1], where_zero_are[0][1:])]
            where_add_text_x = [bins[np.argmax(np.abs(hist_to_plot_fin[am1:am + 1])) + am1] for am1, am in
                                zip(where_zero_are[0][:-1], where_zero_are[0][1:])]
            where_add_text_y = [np.max(np.abs(hist_to_plot_fin[am1:am + 1])) for am1, am in
                                zip(where_zero_are[0][:-1], where_zero_are[0][1:])]
            # integrated_percent_change.append(hist_to_plot_fin[where_zero_are[0][-1]:].sum()/hist_to_plot_fin_ref[where_zero_are[0][-1]:].sum()*100)
        else:
            integrated_percent_change = [
                hist_to_plot_fin[where_zero_are[0][1]:].sum() / hist_to_plot_fin_ref[where_zero_are[0][1]:].sum()]
            where_add_text_x = [bins[np.argmax(np.abs(hist_to_plot_fin[where_zero_are[0][1]:])) + where_zero_are[0][1]]]
            where_add_text_y = [np.max(np.abs(hist_to_plot_fin[where_zero_are[0][1]:]))]
        print integrated_percent_change
        print where_add_text_x
        print where_add_text_y
        for annotx, annoty, value in zip(where_add_text_x, where_add_text_y, integrated_percent_change):
            if (np.isfinite(value)) & (annotx > 0.1) & (np.abs(value) > 6):
                sign = value / np.abs(value)
                add_sign = {-1: '', 1: '+'}
                ax.text(annotx + np.random.rand() * 2, sign * (annoty + np.random.rand() * 0.0002),
                        '{}{:d}\%'.format(add_sign[sign], int(value)), fontsize=10,
                        color=color, )  # ha='center', va='center')

                # ax.annotate('{}{:d}\%'.format(add_sign[sign], int(value)),
                # ax.annotate('{}{:d}\%'.format(add_sign[sign], int(value)),
                #             xy=(annotx, annoty* 2 / 3. * sign), xytext=(annotx, annoty*3/2. * sign),
                #             color=color,
                #             arrowprops=dict(facecolor='black',# shrink=0.005,))
                #                             headwidth=0.5, width=0.1, headlength=0.01))#, width=0.1))
    else:
        integrated_percent_change = percent_change.sum()

    # y_ticks = ax.yaxis.get_major_ticks()
    # y_ticks[5].set_visible(False)
    # y_ticks[6].label.set_visible(False)


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


def mask_3Dcube_with_2Dmask(cube, mask2d):
    '''
    Returns a cube with masked data where the condition cdt is met
    
    :param iris.cube cube:3D cube 
    :param np.ma.maskedarray mask2d: 2D mask to map the 3D cube
    :return iris.cube masked_cube: masked cube
    '''
    masked_cube = cube.copy()
    if not isinstance(masked_cube.data, np.ma.core.MaskedArray):
        masked_cube.data = np.ma.asarray(
            cube.data)  # np.ma.MaskedArray(cube.data,  mask = False)# temperature data below freezing
    mask3d = iris.util.broadcast_to_shape(mask2d, cube.shape, (1, 2))
    masked_cube.data.mask = mask3d | masked_cube.data.mask
    return masked_cube


def generate_year_indices(start_year, end_year, bootstrap_samples_split,
                          bootstrap_no_years, bootstrap_filename, master_ref_year,
                          incomplete_years=False, bootstrap_year_set=[], print_diags=False):

    with open(bootstrap_filename) as bootfile:
        bootfile_content = bootfile.readlines()
        bootfile_content = [x.strip() for x in bootfile_content]
    print('Bootstrapping: read input file')

    all_year_indices = {}

    for loop_boot in np.arange(bootstrap_samples_split[0],
                               bootstrap_samples_split[1]):
        no_years = bootstrap_no_years
        bootstrap_sample_int = np.fromstring(bootfile_content[loop_boot],
                                             dtype=int, count=no_years + 1, sep=' ')
        isample = bootstrap_sample_int[0]
        year_index = bootstrap_sample_int[1:no_years + 1]
        year_index = np.array(year_index)

        first_year = start_year
        last_year = end_year
        if print_diags:
            print('first year: ', first_year)
            print('last year: ', last_year)
        year_label = first_year
        syearref = master_ref_year
        nyearref = bootstrap_no_years

        if incomplete_years:
            print('!!!!! INCOMPLETE YEARS FOR CP4FUT !!!!')
            yrange = bootstrap_year_set
            nyears = len(yrange)
        else:
            nyears = last_year - first_year
            nyears += 1
            yrange = first_year + np.arange(nyears)

        adj_year_index = year_index.copy()

        if incomplete_years:
            if first_year != syearref:
                raise UserWarning('Incomplete years and inconsistent start yr')
            if nyears != nyearref:
                raise UserWarning(
                    'No of years of data ({:d}) is incompatible with bootstrap master ({:d})'.format(nyears, nyearref))
        else:
            if nyears != nyearref:
                raise UserWarning(
                    'No of years of data ({:d}) is incompatible with bootstrap master ({:d})'.format(nyears, nyearref))
            elif first_year > syearref:
                # ensure consistent sampling of overlapping years
                year_offset = first_year - syearref
                adj_year_index[year_index < year_offset] += nyearref
                adj_year_index = adj_year_index - year_offset
                if print_diags:
                    print('Master yr index   ', year_index)
                    print('Adjusted yr index ', adj_year_index)
            elif first_year < syearref:
                # ensure consistent sampling of overlapping years
                year_offset = syearref - first_year
                adj_year_index[year_index > (nyearref - year_offset - 1)] -= nyearref
                adj_year_index = adj_year_index + year_offset
                if print_diags:
                    print('Master yr index   ', year_index)
                    print('Adjusted yr index ', adj_year_index)

        all_year_indices[isample] = [adj_year_index, yrange]

    return all_year_indices


def generate_bootstrapped_hist(dict_by_years, adj_year_index, yrange):

    diag_all = dict_by_years[str(yrange[adj_year_index[0]])].copy()

    for yindex in adj_year_index[1:]:
        year_select = yrange[yindex]
        diag_all = diag_all + dict_by_years[str(year_select)]
    diag_all /= len(adj_year_index)
    return diag_all


def generate_bootstrapped_cubes(dict_by_years, all_year_indices, intensity_bins,
                               duration_bins, following_drought_bins=None):

    cubelist = iris.cube.CubeList()
    for isample in all_year_indices:
        yearlist = all_year_indices[isample]
        data = generate_bootstrapped_hist(dict_by_years,
                                          adj_year_index=yearlist[0],
                                          yrange=yearlist[1])
        cube = hist2cube(data, intensity_bins, duration_bins, following_drought_bins)
        boot_coord = iris.coords.AuxCoord([isample], long_name='bootstrap_sample_number', units='1')
        cube.add_aux_coord(boot_coord)
        cubelist.append(cube)

    bootstrapped_cube = cubelist.merge_cube()

    return bootstrapped_cube


def generate_bootstrapped_event_cubes(dict_by_years, all_year_indices, event, nb_points, season_length):

    cubelist = iris.cube.CubeList()
    for isample in all_year_indices:
        yearlist = all_year_indices[isample]
        data = generate_bootstrapped_hist(dict_by_years,
                                          adj_year_index=yearlist[0],
                                          yrange=yearlist[1])
        _, event_nb = event.freq_or_nb_of_event(data,
                                                nb_points=nb_points,
                                                diff='percdiff', season_length=season_length)
        cube = event2cube(event_nb, )
        boot_coord = iris.coords.AuxCoord([isample], long_name='bootstrap_sample_number', units='1')
        cube.add_aux_coord(boot_coord)
        cubelist.append(cube)

    bootstrapped_cube = cubelist.merge_cube()

    return bootstrapped_cube


### ME2 add this function
def generate_bootstrapped_hist_cubes(dict_by_years, all_year_indices, bins):

    cubelist = iris.cube.CubeList()
    for isample in all_year_indices:
        yearlist = all_year_indices[isample]
        data = generate_bootstrapped_hist(dict_by_years,
                                          adj_year_index=yearlist[0],
                                          yrange=yearlist[1])
        cube = hist2cube(data, bins)
        boot_coord = iris.coords.AuxCoord([isample], long_name='bootstrap_sample_number', units='1')
        cube.add_aux_coord(boot_coord)
        cubelist.append(cube)

    bootstrapped_cube = cubelist.merge_cube()

    return bootstrapped_cube


def generate_cube_from_hist_with_year_dim(dict_by_years, bins):
    cubelist = iris.cube.CubeList()
    for year in dict_by_years:
        cube = hist2cube(dict_by_years[year], bins)
        boot_coord = iris.coords.AuxCoord([int(year)], long_name='years', units='1')#, dtype='|S4')
        cube.add_aux_coord(boot_coord)
        cubelist.append(cube)

    year_cube = cubelist.merge_cube()

    return year_cube


def event2cube(figure):

    data = np.array([figure,])
    event_coord = iris.coords.AuxCoord(1,
                                       long_name='nb_events',
                                       units='1')
    dim_coords_and_dims = ((event_coord, 0), )
    cube = iris.cube.Cube(data, aux_coords_and_dims=dim_coords_and_dims)
    return cube


def hist2cube(data, intensity_bins, duration_bins=None, following_drought_bins=None):

    duration_bin_coord = None
    intens_bin_coord = None
    follow_drougth_coord = None
    if duration_bins is not None:
        duration_bin_coord = iris.coords.AuxCoord(duration_bins,
                                              long_name='duration bins',
                                              units='1')

    if intensity_bins is not None:
        intens_bin_coord = iris.coords.AuxCoord(intensity_bins,
                                              long_name='intensity bins',
                                              units='1')

    if following_drought_bins is not None:
        follow_drougth_coord = iris.coords.AuxCoord(following_drought_bins,
                                                  long_name='following drought bins',
                                                  units='1')

    if duration_bin_coord:
        if intens_bin_coord:
            dim_coords_and_dims = ((intens_bin_coord, 1), (duration_bin_coord, 0))
            if follow_drougth_coord:
                dim_coords_and_dims = ((intens_bin_coord, 0), (duration_bin_coord, 1), (follow_drougth_coord, 2))
        else:
            dim_coords_and_dims = ((duration_bin_coord, 0), )
    else:
        dim_coords_and_dims = ((intens_bin_coord, 0), )

    cube = iris.cube.Cube(data, aux_coords_and_dims=dim_coords_and_dims)#, aux_coords_and_dims=aux_coords_and_dims)
    return cube

def estimate_local_pvalue_from_bootstrap(control_cube, control_cube_simulated, future_cube,
                                         future_cube_simulated, f_ts=iam.subtract, print_statements=False):
    '''
    returns p_val as the proportion of sample differences
    (future_cube_simulated-control_cube_simulated)
    lying outside the difference between future-control.

    :param iris.cube control_cube: 2D cube
    :param iris.cube control_cube_simulated: 3D cube with 'simulated_sample' coordinate
    :param iris.cube future_cube: 2D cube
    :param iris.cube future_cube_simulated: 3D cube with 'simulated_sample' coordinate
    :param func f_ts: default: iris.analysis.math.subtract

    :return p_val: cube with pvalues
    '''

    ''' Do not standardize ts because the distribution of ts_simulated is not standardized? '''
    ts = f_ts(control_cube, future_cube)  # n dimension
    resampled_ts = control_cube_simulated.copy()
    resampled_ts.data = future_cube_simulated.data - control_cube_simulated.data
    resampled_ts_mean = resampled_ts.collapsed('bootstrap_sample_number', iana.MEAN)
    null_ts = resampled_ts - resampled_ts_mean  # n+1 dimension
    min_pval = 1.0 / len(resampled_ts.coord('bootstrap_sample_number').points)
    if print_statements:
        print(min_pval)
        print(resampled_ts.data.shape, ts.data.shape)
    _, ts_ = np.broadcast_arrays(resampled_ts.data, ts.data[None, ...])

    p_ = null_ts.copy()
    p_.data = (ts_ > null_ts.data).astype(float)  # Slooooooooooooooooow
    # print('Arghhhhh. The wall is collapsing!')

    p_val = p_.collapsed('bootstrap_sample_number', iana.COUNT,  # iana.PROPORTION doesn't always work...
                         function=lambda values: values > 0.5)
    if print_statements:
        print('p_val 000')
        print(p_val.data.shape, p_val.data)

    p_val = p_val / np.float(len(resampled_ts.coord('bootstrap_sample_number').points))

    if print_statements:
        print('p_val orig:')
        print(p_val.data.max(), p_val.data.min(), np.median(p_val.data), \
          np.percentile(p_val.data, [1.0, 5.0, 10.0, 25.0, 75.0, 90.0, 95.0, 99.0]))

    #### two tailed test ######
    p_val.data[p_val.data > 0.5] = 1.0 - p_val.data[p_val.data > 0.5]

    if print_statements:
        print('p_val two tailed:')
        print(p_val.data.max(), p_val.data.min(), np.median(p_val.data), \
              np.percentile(p_val.data, [1.0, 5.0, 10.0, 25.0, 75.0, 90.0, 95.0, 99.0]))

    p_val.data[p_val.data <= min_pval] = min_pval

    if print_statements:
        print('p_val final:')
        print(p_val.data.max(), p_val.data.min(), np.median(p_val.data), \
              np.percentile(p_val.data, [1.0, 5.0, 10.0, 25.0, 75.0, 90.0, 95.0, 99.0]))

    return p_val.data

def add_CORDEX_grid(cube):
    if 'regridded_to' in cube.attributes:
       nc = netCDF4.Dataset('/home/users/sberthou/PrecipExtremes/{}_grid.nc'.format(cube.attributes['regridded_to']), 'r')
    else:
       nc = netCDF4.Dataset('/home/users/sberthou/PrecipExtremes/EUROCORDEX_grid.nc', 'r')
       # raise UserWarning('No info found for grid, make sure that the history attributes contains cdo remapcon and CORDEX')
    lat_tab = nc.variables['rlat'] 
    lon_tab = nc.variables['rlon']
    pole = nc.variables['rotated_latitude_longitude'] 
    coord_sys = iris.coord_systems.RotatedGeogCS(pole.grid_north_pole_latitude, 
                                                 pole.grid_north_pole_longitude, 
                                                 ellipsoid=iris.coord_systems.GeogCS(6371229.0))
    lat_coord = iris.coords.DimCoord(lat_tab, standard_name = 'grid_latitude', 
                                              units=cf_units.Unit('degrees'), 
                                              coord_system = coord_sys)
    lon_coord = iris.coords.DimCoord(lon_tab, standard_name = 'grid_longitude', 
                                              units=cf_units.Unit('degrees'), 
                                              coord_system = coord_sys)
    if len(cube.shape) == 3:
        cube.add_dim_coord(lat_coord, 1)
        cube.add_dim_coord(lon_coord, 2)
    else:
        cube.add_dim_coord(lat_coord, 0)
        cube.add_dim_coord(lon_coord, 1)


def handle_ALADIN_Lambert_grid(cube):
   coord_system = cube.coord_system()
   assert isinstance(coord_system, LambertConformal)
   
   keys = ('central_lat', 'central_lon', 'false_easting', 'false_northing', 'secant_latitudes')
   kwargs = {key:getattr(coord_system, key) for key in keys}
   
   # Replace the secant latitudes 1-tuple with a 2-tuple == two of the same value.
   secant_lats = kwargs['secant_latitudes']
   assert len(secant_lats) == 1
   secant_lat = secant_lats[0]
   kwargs['secant_latitudes'] = (secant_lat, secant_lat)
   kwargs['false_easting'] = 2827740    
   kwargs['false_northing'] = 2816550
   # Create a new coord system from the modified keys.
   new_coord_system = LambertConformal(**kwargs)
   
   # Replace the coord system of the grid coordinates.
   # ALSO: convert units to 'm', which seems to be needed for Cartopy plotting ?
   for dim_name in ('projection_x_coordinate', 'projection_y_coordinate'):
      coord = cube.coord(dim_name)
      coord.coord_system = new_coord_system
      # Also fix projection coordinate units : Seems to be wrong in the file ??
      assert coord.units == 'km'
      coord.convert_units('m')
      assert coord.units == 'm'


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


####update11/10
def plot_map(subregion_files, bigregion, srex_reg=None, reg_todo = None):
    fig = plt.figure()
    subpl = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(), frameon=False)
    # subpl.stock_img()
    subpl.coastlines('50m')
    # if not reg_todo:
    subpl.add_feature(cfeature.BORDERS, linewidth=0.5)
    # gl = subpl.gridlines(draw_labels=True)
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    for ax_pos in ['top', 'bottom', 'left', 'right']:
        subpl.spines[ax_pos].set_visible(False)
    if not reg_todo:
        # pass
        print subregion_files
        col_order = plot_masked_subregions(subregion_files, subpl)
    else:
        for r in reg_todo:
            x, y = srex_reg[r].exterior.xy
            subpl.plot(x, y, color='k')  # plot regional box lines
            xc, yc = np.array(srex_reg[r].exterior.centroid, )
            if r == 'CSA':
                r = 'SSZ'
            subpl.text(xc, yc, r, horizontalalignment='center',
                       verticalalignment='center', color='red',
                       fontdict={'weight': 'bold'})  # add region codes to boxes
    subpl.set_extent(bigregion, ccrs.PlateCarree())
    return subpl


def plot_pie_chart(ax, percent_sig_per_season, reg, comp, comp_let=[], percent_interval=0.5):
    size = 0.3
    vals = [1, 1, 1, 1, ]
    thresholds_col = (0, percent_interval, 1.0)#(0, 0.5, 0.9, 1.0)
    total_colors = len(thresholds_col) - 1
    cmaps = {'jja': plt.get_cmap("Reds", total_colors+1),
             'djf': plt.get_cmap("Blues", total_colors+1),
             'son': plt.get_cmap("Oranges", total_colors+1),
             'mam': plt.get_cmap("Greens", total_colors+1)}

    intensities = ['60to400', '10to60', '1to10', ]
    seasons = ['jja', 'mam', 'djf', 'son']

    for i, intensity in enumerate(intensities):
        outer_colors = []
        letters = []
        for season in seasons:
            num_col = percent_sig_per_season[season][reg][comp][intensity]
            colcol = total_colors - 1 
            for ii in range(total_colors-1, 0, -1): 
                if num_col < thresholds_col[ii]: 
                    colcol = ii-1 
            
            color = cmaps[season]([colcol,])[0]
            outer_colors.append(color)
            num_let = ''
            if comp_let:
               if colcol > 0:
                  print 'assuming comp_let[0] is CORDEX, comp_let[1] is PRIMAVERA'
                  cordex_dist = percent_sig_per_season[season][reg][comp_let[0]][intensity]
                  prim_dist = percent_sig_per_season[season][reg][comp_let[1]][intensity]
                  if prim_dist < cordex_dist - 0.1:
                     num_let = 'P'
                  elif prim_dist > cordex_dist + 0.1 :
                     num_let = 'C'
               else:
                  if percent_sig_per_season[season][reg][comp_let[0]][intensity] < 0.7:
                     if percent_sig_per_season[season][reg][comp_let[1]][intensity] < 0.7:
                        num_let = '='
            letters.append(num_let)           
        print letters 
        _, texts = ax.pie(vals, radius=1-i*size, colors=outer_colors, labels=letters, labeldistance=0.75,
               wedgeprops=dict(width=size, ), startangle=45)
        ax.set_aspect('equal')
        for tt in texts:
           tt.set_fontsize(5)
           tt.set_color('w')
           tt.set_fontweight('bold')

def plot_pie_inset(subregion_files, ax, width, percent_sig_per_season, figdir, comp_col='CORDEX50vsPRIMAVERA', comp_let=[], reg_todo=None, plot_type='', distrib='exponential', percent_interval=0.5, sig_val=1):
    if not reg_todo:
        european_centroids = {'IP': (-4, 40),
                              'FR': (1, 46),
                              'BI': (-2, 53),
                              'CE': (11, 53),
                              'SC': (12, 62),
                              'AL': (9, 46),
                              'MD': (14, 40),
                              'NEE': (24, 54),
                              'CA': (22, 47),
                              }
        for reg in subregion_files:
            ilon, ilat = european_centroids[reg]

            ax_sub = inset_axes(ax, width=width, height=width, loc=10,
                                bbox_to_anchor=(ilon, ilat),
                                bbox_transform=ax.transData,
                                borderpad=0)
            plot_pie_chart(ax_sub, percent_sig_per_season, reg, comp_col, comp_let, percent_interval)
    savename = '{}_{}_pie_chart_{}_{}_{}_{}_{}'.format('_'.join([r for r in subregion_files]), comp_col, distrib, plot_type, percent_interval, sig_val, int(np.random.random()*100))
    plt.savefig(figdir+'{}.png'.format(savename), dpi=300)
    plt.savefig(figdir+'{}.pdf'.format(savename))

# def generate_global_pvalue_mask(control_cube, control_cube_simulated, future_cube,
#                                 future_cube_simulated, sig_lev=0.05):
#     '''
#     This function calculates the local pvalues from the bootstrapping
#     (calling estimate_local_pvalue_from_bootstrap), and then calculates the global p_value with
#     multicomp.multipletests and generates a mask of significant values
#     (then can be applied to the diff cube and plotted using hatches = ["..."] )
#
#     :param iris.cube control_cube: 2D cube
#     :param iris.cube control_cube_simulated: 3D cube with 'simulated_sample' coordinate
#     :param iris.cube future_cube: 2D cube
#     :param iris.cube future_cube_simulated: 3D cube with 'simulated_sample' coordinate
#     :param float sig_lev: global level of significance wanted to be achieved
#
#     :return annual_sig_test_mask
#     '''
#     alpha_fdr_adj = 2.0
#
#     p_val = estimate_local_pvalue_from_bootstrap(control_cube, control_cube_simulated,
#                                                  future_cube, future_cube_simulated)
#     p_val_compressed = p_val.data.compressed() if isinstance(p_val.data, np.ma.MaskedArray) else p_val.data.flatten()
#     p_val2 = multicomp.multipletests(p_val_compressed,
#                                      alpha=alpha_fdr_adj * sig_lev / 2.0,
#                                      method='fdr_bh', is_sorted=False,
#                                      returnsorted=False)
#     annual_sig_test_mask_raw = np.logical_not(p_val2[0])
#
#     if isinstance(p_val.data, np.ma.MaskedArray):  # reshape and remask
#         annual_sig_test_mask = np.empty_like(p_val.data)
#         np.place(annual_sig_test_mask, annual_sig_test_mask.mask, 1)
#         np.place(annual_sig_test_mask, ~annual_sig_test_mask.mask,
#                  annual_sig_test_mask_raw)
#     else:
#         annual_sig_test_mask = np.reshape(annual_sig_test_mask_raw, p_val.data.shape)
#
#     print(p_val2[0][0:10], p_val2[1][0:10], p_val_compressed[0:10])
#     print(p_val.data.max(), p_val.data.min(), p_val.data.shape)
#     print(p_val2[1].max(), p_val2[1].min(), np.median(p_val2[1]), np.mean(p_val2[1]), p_val2[1].shape)
#
#     return annual_sig_test_mask
