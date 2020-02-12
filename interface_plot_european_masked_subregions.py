'''
Created on Jun 29, 2016

@author: cbates

For plotting the interface and the histogram plots to be visible
'''

import __init__
import european_masked_subregion as ems
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", type=str, help="for which country? (Alps/France/Spain/ukcpobs/nimrod/mali/niger/benin")
parser.add_argument("-r", type=str, help="for the files regridded on n512 (or model for ammacatch ", default='False')
parser.add_argument("-f", type=str, help="", default='d')
parser.add_argument("-o", type=str, help="Other (ukcpm_runs, ammacatch", default='')
parser.add_argument("-s", nargs='*', type=str, help="Season list", default=[None])
parser.add_argument("-b", type=str, help="traditional or/and martin distribution?", default=['traditional'])
parser.add_argument("-n", type=str, help="normed / frac / normed wet / contrib", default='normed')
parser.add_argument("-g", type=str, help="1d or 2d?", default='1d')
parser.add_argument("-w", type=str, help="wet or dry?", default='wet')
parser.add_argument("-t", type=float, help="threshold (mm/day?)?", default=1.0)

args = parser.parse_args()

if args.c is None or args.r is None or args.f is None or args.o is None or args.s is None or args.b is None:
    parser.print_help()
    print "\nPlease provide valid arguments\n"
    exit()


def bool_from_string(stringb):
    if stringb.lower() == "true":
        bools = True
    elif stringb.lower() == "false":
        bools = False
    else:
        raise ValueError('string must be "True" or "False"')
    return bools

country = args.c
n512 = bool_from_string(args.r)
frequency = args.f
other = args.o
seasons = args.s
bintype = args.b
fractional = args.n
distribution = args.b
hist_type = args.g
wet_or_dry_spells = args.w
thrs = args.t

### ME2: this defines ensembles of models: if a model or dataset is not in an ensemble,
# the lower/higher bounds of the bootstrapped interval will come from inter-annual variability in the pdf.
# e.g. for the observations.
ensemble_dict = {'PRIMAVERA': ['MPI-ESM1-2-XR_r1i1p1f1_EUROCORDEX',
                               'HadGEM3-GC31-HM_r1i1p1f1_EUROCORDEX',
		 	       'CNRM-CM6-1-HR_r1i1p1f2_EUROCORDEX',
		 	       'EC-Earth3P-HR_r1i1p2f1_EUROCORDEX',
			       'CMCC-CM2-VHR4_r1i1p1f1_EUROCORDEX', 
		 	       'ECMWF-IFS-HR_r1i1p1f1_EUROCORDEX',
                               #'ECMWF-IFS-HR_r2i1p1f1_EUROCORDEX',
                               ###'ECMWF-IFS-HR_r3i1p1f1_EUROCORDEX',
                               #'ECMWF-IFS-HR_r4i1p1f1_EUROCORDEX',
                               #'ECMWF-IFS-HR_r5i1p1f1_EUROCORDEX',
                               #'ECMWF-IFS-HR_r6i1p1f1_EUROCORDEX',
		 	       ],
                 #'PRIMAVERA': ['MPI-ESM1-2-XR_r1i1p1f1',
                 #              'HadGEM3-GC31-HM_r1i1p1f1',
	         #	       'CNRM-CM6-1-HR_r1i1p1f2',
		 #              'EC-Earth3P-HR_r1i1p2f1',
		 #	       'CMCC-CM2-VHR4_r1i1p1f1',
	         #               'ECMWF-IFS-HR_r1i1p1f1',
                 #               'ECMWF-IFS-HR_r2i1p1f1',
                 #               'ECMWF-IFS-HR_r3i1p1f1',
                 #               'ECMWF-IFS-HR_r4i1p1f1',
                 #               'ECMWF-IFS-HR_r5i1p1f1',
                 #               'ECMWF-IFS-HR_r6i1p1f1',
		 #               ],
		 #'PRIMAVERA LR': ['EC-Earth3P_r1i1p2f1',
		 #		  'MPI-ESM1-2-HR_r1i1p1f1',
		 #		  'CNRM-CM6-1_r1i1p1f2',
		 #		  'HadGEM3-GC31-LL_r1i1p1f1',
		 #		  'HadGEM3-GC31-LL_r1i2p1f1',
		 #		  'HadGEM3-GC31-LL_r1i3p1f1',
		 #		  'HadGEM3-GC31-LL_r1i4p1f1',
		 #		  'HadGEM3-GC31-LL_r1i5p1f1',
		 #		  'HadGEM3-GC31-LL_r1i6p1f1',
		 #		  'HadGEM3-GC31-LL_r1i7p1f1',
		 #		  'HadGEM3-GC31-LL_r1i8p1f1',
		 #                'CMCC-CM2-HR4_r1i1p1f1',
		 # 		  'ECMWF-IFS-LR_r1i1p1f1',
		 #		  'ECMWF-IFS-LR_r2i1p1f1',
		 #		  'ECMWF-IFS-LR_r3i1p1f1',
		 #		  'ECMWF-IFS-LR_r4i1p1f1',
		 #		  'ECMWF-IFS-LR_r5i1p1f1',
		 #		  'ECMWF-IFS-LR_r6i1p1f1',
		 #		  'ECMWF-IFS-LR_r7i1p1f1',
		 #		  'ECMWF-IFS-LR_r8i1p1f1',
		 # 		  ],
		 #'CORDEX12': ['CCLM-MPI',
                 #	   'CCLM-CNRM',
                 #	   'CCLM-EC-EARTH',
                 #	   'HIRHAM-CNRM',  
                 #	   'HIRHAM-EC-EARTH',		     
                 #	   'RACMO-CNRM',   
                 #	   'RACMO-EC-EARTH',		     
                 #	   'RACMO-HadGEM', 
                 #	   'REMO-MPI',
                 #	   'RCA-CNRM',
                 #	   'RCA-EC-EARTH', 
                 #	   'RCA-HadGEM',   
                 #	   'RCA-MPI',
                 #	   #'HIRHAM-NorESM', #select only RCM forced by CMIP5 GCMs similar to PRIMAVERA
                 #	   #'REMO-IPSL',
                 #	   #'REMO-NorESM',  
                 #	   #'REMO-GFDL',
                 #	   #'WRF-IPSL',  #zero/nan values in json file
                 #	   #'RACMO-NorESM', 
                 #	   #'RCA-IPSL',
                 #	   #'RCA-NorESM',
		 #	   ],
		 'CORDEX-44':['CCLM-MPI_50km',
                 	   'CCLM5-CNRM_50km',		      
                 	   'CCLM5-EC-EARTH_50km',	     
                 	   'CCLM5-MPI_50km',
		 	           'CCLM5-HadGEM_50km', 
                 	   'HIRHAM-EC-EARTH_50km',	   
                 	   'WRF-EC-EARTH_50km', 
                 	   'RACMO-EC-EARTH_50km',	   
                	   'RACMO-HadGEM_50km', 	   
                 	   'REMO-MPI_50km',
                 	   'RCA-CNRM_50km',
                 	   'RCA-EC-EARTH_50km', 	   
                 	   'RCA-HadGEM_50km',		   
                 	   'RCA-MPI_50km',
                 	   'CCLM5-MIROC_50km',  #only include RCMs forced with CMIP5 GCMs similar to PRIMAVERA  	   
                 	   'WRF-IPSL_50km',  #zero/nan values in json file: problems when creating json
                 	   'RCA-CanESM_50km', #		  
                 	   'RCA-CSIRO_50km',	#	  
                 	   'RCA-IPSL_50km',#
                 	   'RCA-MIROC_50km',#		  
                 	   'RCA-NorESM_50km',#
                 	   'RCA-GFDL_50km',#
                       'RegCM-HadGEM_50km',
                       'ALADIN-CNRM1_50km',
                       'ALARO-CNRM_50km',
                       'WRF-CanESM_50km',#
                       #'REMO-MPI2_50km', # not including second REMO member
		 	  ],
	         #'CMIP5':['MPI-ESM-LR',
                 #	  'CNRM-CM5',
                 #	  'EC-EARTH',
                 #	  'HadGEM2-ES',
                 #	  #'NorESM1-M', #select only CMIP5 GCMs similar to PRIMAVERA
                 #	  #'IPSL-CM5A-LR',
                 #	  #'IPSL-CM5A-MR',
                 #	  #'GFDL-ESM2G',
                 #	  #'MIROC5',
                 #	  #'CanESM2',
                 #	  #'CSIRO-Mk3-6-0',
                 #	  #'GFDL-ESM2M',
		 #	  ],
		 }
#what_to_bootstrap = ['PRIMAVERA', 'CORDEX-44', 'obs_on_12km']
what_to_bootstrap = ['PRIMAVERA', 'CORDEX-44', 'obs_cordex50']
#what_to_bootstrap = [ 'PRIMAVERA']


bootstrap_config = {
                   'PRIMAVERA':
			{'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
			 'bootstrap_no_years': len(ensemble_dict['PRIMAVERA']),  # number of PRIMAVERA models
			 'bootstrap_sample_size': 99,  # NB if this goes > 999 (3 digits) formatting fails
 		 'bootstrap_samples_split': [0, 99],  # NB if this goes > 999 (3 digits) formatting fails
			 'bootstrap_filename': 'bootstrap_ranint_master_primavera_100.txt',  # I/O name for bootstrap years
			 'master_ref_year': 1997, # random
			},
		    #'PRIMAVERA LR':
		    #	{'bootstrap_produce_file': True,  # True if want to generate bootstrapping random years rather than use exisiting file
		    #	 'bootstrap_no_years': len(ensemble_dict['PRIMAVERA LR']),  # number of PRIMAVERA models
		    #	 'bootstrap_sample_size': 999,  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_samples_split': [0, 999],  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_filename': 'bootstrap_ranint_master_primavera_lr.txt',  # I/O name for bootstrap years
		    #	 'master_ref_year': 1997, # random
		    #	},
		    #'CMIP5':
		    #	{'bootstrap_produce_file': True,  # True if want to generate bootstrapping random years rather than use exisiting file
		    #	 'bootstrap_no_years': len(ensemble_dict['CMIP5']),  # number of PRIMAVERA models
		    #	 'bootstrap_sample_size': 999,  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_samples_split': [0, 999],  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_filename': 'bootstrap_ranint_master_cmip5.txt',  # I/O name for bootstrap years
		    #	 'master_ref_year': 1997, # random
		    #	},
		    #'CORDEX12':
		    #	{'bootstrap_produce_file': True,  # True if want to generate bootstrapping random years rather than use exisiting file
		    #	 'bootstrap_no_years': len(ensemble_dict['CORDEX12']),  # number of PRIMAVERA models
		    #	 'bootstrap_sample_size': 999,  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_samples_split': [0, 999],  # NB if this goes > 999 (3 digits) formatting fails
		    #	 'bootstrap_filename': 'bootstrap_ranint_master_cordex12.txt',  # I/O name for bootstrap years
		    #	 'master_ref_year': 1997, # random
		    #	},
		    'CORDEX-44':
		       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		        'bootstrap_no_years': len(ensemble_dict['CORDEX-44']),  # number of PRIMAVERA models
		        'bootstrap_sample_size': 999,  # NB if this goes > 999 (3 digits) formatting fails
		        'bootstrap_samples_split': [0, 999],  # NB if this goes > 999 (3 digits) formatting fails
		        'bootstrap_filename': 'bootstrap_ranint_master_cordex50.txt',  # I/O name for bootstrap years
		        'master_ref_year': 1997, # random
		       },		       
		    #'obs_on_12km':
		    'obs_cordex50':
		       {
		       'run_IP_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
		       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		        'bootstrap_no_years': 33,  # nb of OBS years, different by region
		        'bootstrap_sample_size': 999,
		        'bootstrap_samples_split': [0, 999],
		        'bootstrap_filename': 'bootstrap_ranint_master_ibp.txt',  # change file name by region
		        'master_ref_year': 1971, # start year, different by region
		       },
                       'run_AL_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_alps.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },		      
                       'run_BI_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999, 
                        'bootstrap_samples_split': [0, 999], 
                        'bootstrap_filename': 'bootstrap_ranint_master_uk.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       }, 
		       'run_CA_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_carpathian.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
		       'run_FR_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                    	'bootstrap_sample_size': 999,
                    	'bootstrap_samples_split': [0, 999],
                    	'bootstrap_filename': 'bootstrap_ranint_master_france.txt',  # change file name by region
                    	'master_ref_year': 1971, # start year, different by region
                       },
                       'run_CE_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_germany.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_MD_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_mediterranean.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_SC_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_scandinavia.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_NEE_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': 999,
                        'bootstrap_samples_split': [0, 999],
                        'bootstrap_filename': 'bootstrap_ranint_master_scandinavia.txt',  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       }, 
		    }   
#ME2 add this:
if bootstrap_config:
    for ens in what_to_bootstrap:
	#if ens == 'obs_on_12km':
	print 'ens = ', ens
	if ens == 'obs_cordex50':
	    for region in bootstrap_config[ens].keys():
	        #print 'region = ', region
		if bootstrap_config[ens][region]['bootstrap_produce_file']:
		    print bootstrap_config[ens][region]['bootstrap_produce_file']
		    with open(bootstrap_config[ens][region]['bootstrap_filename'], "w") as file:
			ran_out = np.random.randint(0, bootstrap_config[ens][region]['bootstrap_no_years'],
						    size=(bootstrap_config[ens][region]['bootstrap_sample_size'],
							  bootstrap_config[ens][region]['bootstrap_no_years']))
			multi_placeholder = "{:2d} " * bootstrap_config[ens][region]['bootstrap_no_years']
			string_placeholder = "{:3d} " + multi_placeholder + "\n"
			for a in range(0, bootstrap_config[ens][region]['bootstrap_sample_size']):
			    file.write(string_placeholder.format(a, *ran_out[a]))
	else:
	    if bootstrap_config[ens]['bootstrap_produce_file']:
		with open(bootstrap_config[ens]['bootstrap_filename'], "w") as file:
		    ran_out = np.random.randint(0, bootstrap_config[ens]['bootstrap_no_years'],
						size=(bootstrap_config[ens]['bootstrap_sample_size'],
						      bootstrap_config[ens]['bootstrap_no_years']))
		    multi_placeholder = "{:2d} " * bootstrap_config[ens]['bootstrap_no_years']
		    string_placeholder = "{:3d} " + multi_placeholder + "\n"
		    for a in range(0, bootstrap_config[ens]['bootstrap_sample_size']):
			file.write(string_placeholder.format(a, *ran_out[a]))

xlim = {'d': (0.9, 500),
        'h': (0.02, 60)}  # (1./24./4., 100)
savefig = True
localdir = '/home/users/sberthou/PrecipExtremes/'  ### ME: must be the same as in process_data_europe, config['datadir']
dict_freq = {
    'h': 24.,
    'd': 1.,
    '3h': 8.,
    '15mn': 24. * 4.
}

json_freq = {'d': '/',
             '3h': '_3hourly/',
             'h': '_hourly/',
             '15mn': '_15mn/'}


class Region(object):
    '''
    This just acts as a dictionary to store many arguments: 
    '''

    def __init__(self, name, lon_start, lon_end, lat_start, lat_end, dictppt):
        self.name = name  # str
        self.lon_start = lon_start  # float
        self.lon_end = lon_end  # float
        self.lat_start = lat_start  # float
        self.lat_end = lat_end  # float
        self.plotppt = dictppt  # dictionary with plotting properties, used for figure = True (see Europe example)


class Alps(object):
    def __init__(self):
        self.region = Region('ALPS', lon_start=3.0, lon_end=18.0, \
                             lat_start=42, lat_end=50, \
                             dictppt={})


class Spain(object):
    def __init__(self):
        self.region = Region('SPAIN', lon_start=-10.0, lon_end=5.0, \
                             lat_start=35.0, lat_end=44.0, \
                             dictppt={})


class France(object):
    def __init__(self):
        self.region = Region('FRANCE', lon_start=-6.0, lon_end=10.0, \
                             lat_start=40, lat_end=52, \
                             dictppt={})


class UK(object):
    def __init__(self):
        self.region = Region('UK', lon_start=-11.0, lon_end=3.0, \
                             lat_start=48, lat_end=62, \
                             dictppt={})


class Germany(object):
    def __init__(self):
        self.region = Region('Germany', lon_start=5.5, lon_end=15.5, \
                             lat_start=47, lat_end=55, \
                             dictppt={})


class Swiss(object):
    def __init__(self):
        self.region = Region('Switzerland', lon_start=5.5, lon_end=11.0, \
                             lat_start=45, lat_end=48, \
                             dictppt={})


class Europe(object):
    def __init__(self):
        self.region = Region('EUROPE',
                             lon_start=-13.0, lon_end=32.0,
                             lat_start=33.0, lat_end=71.0, #63
                             dictppt={})


class Niger(object):
    def __init__(self):
        self.region = Region('NIGER', lon_start=-6, lon_end=6.0, \
                             lat_start=4.5, lat_end=20, \
                             dictppt={})


_, maskedreglist = ems.get_runlist_region(frequency, n512, country, other)

country_limits = {'spain': Spain().region,
                  'ukcpobs': UK().region,
                  'nimrod': UK().region,
                  'france': France().region,
                  'switzerland': Swiss().region,
                  'germany': Germany().region,
                  'alps': Alps().region,
                  'niger': Niger().region,
                  'benin': Niger().region,
                  'mali': Niger().region,
                  'med': Europe().region,
                  'cen_europe': Europe().region,
                  'above_1500': Europe().region,
                  'med': Europe().region,
                  'neur': Europe().region,
                  'europe': Europe().region,
                  'mediterranean': Europe().region,
                  'prudence': Europe().region,
                  }

region_limits = country_limits[country]
bigregion = (region_limits.lon_start, region_limits.lon_end,
             region_limits.lat_start, region_limits.lat_end)

subregion_files = maskedreglist
print maskedreglist
if hist_type == '1d':
    hists = ['run_{}_{}.json'.format(reg_masked, distribution) for reg_masked in maskedreglist.keys()]
else:
    hists = ['hist2d_{}_{}_{:d}.json'.format(reg_masked, wet_or_dry_spells, int(thrs)) for reg_masked in
             maskedreglist.keys()]

if other == 'ukcpm_runs':
    config = {'jsondir': 'srex_pdf_json/',
              'histfile': {country: hists,  # country
                           '2p2erai_on_nimrod': hists,
                           # '2p2erai_on_ukcpobs':hists,
                           #        '2p2':hists,
                           'UKV_4km': hists,
                           'UKV_2.2km_14': hists,
                           # 'UKV_2p2km_24':hists,
                           'UKV_2p2km_var': hists,
                           'UKV_1p5km': hists,
                           'UKV_4km_nest12': hists,
                           }
              }

elif other == 'ammacatch':
    if frequency == '3h':
        config = {'jsondir': 'srex_pdf_json' + json_freq[frequency],
                  'histfile': {
                      'CP4': hists,
                      'amma-catch': hists,
                      'TRMM': hists,
                      'CMORPH': hists,
                  }
                  }
    else:
        print 'amma-catch h'
        config = {'jsondir': 'srex_pdf_json_hourly_jun17/',  # +json_freq[frequency],
                  'histfile': {
                      # 'cp4':hists,
                      'amma-catch': hists,
                  }
                  }


elif n512:
    config = {'jsondir': 'srex_pdf_json' + json_freq[frequency],
              'histfile': {  # country+'_on_n512':hists,
                  '2p2_on_n512': hists,
                  'n512': hists,
                  'n512_future': hists,
                  '2p2f_on_n512': hists,
              }
              }

else:

    config = {
              'jsondir': '{}/srex_pdf_json{}'.format(localdir, json_freq[frequency]),
              'histfile': {# MED: this is where you need to put CORDEX model names for plotting (name must be the same as the one in european_masked_subregion.py)
                           #
			   #EUR11 models
			   #'CCLM-MPI': hists,
                           #'CCLM-CNRM': hists,
                           #'CCLM-EC-EARTH': hists,
                           #'HIRHAM-CNRM': hists,  
                           #'HIRHAM-EC-EARTH': hists,		       
                           #'RACMO-CNRM': hists,   
                           #'RACMO-EC-EARTH': hists,		       
                           #'RACMO-HadGEM': hists, 
                           #'REMO-MPI': hists,
                           #'RCA-CNRM': hists,
                           #'RCA-EC-EARTH': hists, 
                           #'RCA-HadGEM': hists,   
                           #'RCA-MPI': hists,
                           #'RCA-NorESM': hists,  #only include RCMs forced with CMIP5 GCMs similar to PRIMAVERA	 
                           #'HIRHAM-NorESM': hists,
                           #'REMO-IPSL': hists,
                           #'REMO-NorESM': hists,  
                           #'REMO-GFDL': hists,
                           #'WRF-IPSL': hists,
                           #'RACMO-NorESM': hists, 
                           #'RCA-IPSL': hists,
                           ###'CCLM-HadGEM': hists, #These have no lon/lat coord  
                           ###'ALADIN-CNRM1': hists,
                           ###'ALADIN-CNRM2': hists,
                           ###'ALARO-CNRM': hists,  
                           #
                           ##EUR44 models
                           'CCLM-MPI_50km': hists,
                           'CCLM5-CNRM_50km': hists,		     
                           'CCLM5-EC-EARTH_50km': hists,	     
                           'CCLM5-MPI_50km': hists,		    
                           'HIRHAM-EC-EARTH_50km': hists,	    
                           'WRF-EC-EARTH_50km': hists, #zero/nan values in json file
                           'RACMO-EC-EARTH_50km': hists,	    
                           'RACMO-HadGEM_50km': hists,  	    
                           'REMO-MPI_50km': hists,
                           'RCA-CNRM_50km': hists,
                           'RCA-EC-EARTH_50km': hists,  	    
                           'RCA-HadGEM_50km': hists,		    
                           'RCA-MPI_50km': hists, 
                           'CCLM5-MIROC_50km': hists,  #only include RCMs forced with CMIP5 GCMs similar to PRIMAVERA	    
                           'WRF-IPSL_50km': hists, #zero/nan values in json file
                           'RCA-CanESM_50km': hists,#		   
                           'RCA-CSIRO_50km': hists,#		   
                           'RCA-IPSL_50km': hists,#
                           'RCA-MIROC_50km': hists,#		   
                           'RCA-NorESM_50km': hists,#
                           'RCA-GFDL_50km': hists,#
                           'CCLM5-HadGEM_50km': hists, #These have no lon/lat coord  	       
                           'ALADIN-CNRM1_50km': hists,	      
                           ###'ALADIN-CNRM2_50km': hists,#		 
                           'WRF-CanESM_50km': hists, 	#	 
                           'RegCM-HadGEM_50km': hists,		 
                           'ALARO-CNRM_50km': hists, 		 
                           #
                           #CMIP5 models
                           #'MPI-ESM-LR': hists,
                           #'CNRM-CM5': hists,
                           #'EC-EARTH': hists,
                           #'HadGEM2-ES': hists,
                           #'NorESM1-M': hists,
                           #'IPSL-CM5A-LR': hists,
                           #'IPSL-CM5A-MR': hists,
                           #'GFDL-ESM2G': hists,
                           #'MIROC5': hists,
                           #'CanESM2': hists,
                           #'CSIRO-Mk3-6-0': hists,
                           #'GFDL-ESM2M': hists,
			   #
			   #HighResMIP high-resolution models (HR)
			   #'EC-Earth3P-HR_r1i1p2f1': hists,   #36km
			   #'HadGEM3-GC31-HM_r1i1p1f1': hists, #25km
			   #'CMCC-CM2-VHR4_r1i1p1f1': hists,   #18km
			   #'MPI-ESM1-2-XR_r1i1p1f1': hists,   #34km
			   #'CNRM-CM6-1-HR_r1i1p1f2': hists,   #50km
			   #'ECMWF-IFS-HR_r1i1p1f1': hists,    #25km regridded at 50km
			   #'ECMWF-IFS-HR_r2i1p1f1': hists,   #put only 1 member per model
			   ##'ECMWF-IFS-HR_r3i1p1f1': hists,
			   #'ECMWF-IFS-HR_r4i1p1f1': hists,
			   #'ECMWF-IFS-HR_r5i1p1f1': hists,
			   #'ECMWF-IFS-HR_r6i1p1f1': hists,
			   #HighResMIP high-resolution models (HR) regridded at 50km
			   'EC-Earth3P-HR_r1i1p2f1_EUROCORDEX': hists, 
			   'HadGEM3-GC31-HM_r1i1p1f1_EUROCORDEX': hists,
			   'CMCC-CM2-VHR4_r1i1p1f1_EUROCORDEX': hists,
			   'MPI-ESM1-2-XR_r1i1p1f1_EUROCORDEX': hists,
			   'CNRM-CM6-1-HR_r1i1p1f2_EUROCORDEX': hists,
			   'ECMWF-IFS-HR_r1i1p1f1_EUROCORDEX': hists,
			   #'ECMWF-IFS-HR_r2i1p1f1_EUROCORDEX': hists,
			   ###'ECMWF-IFS-HR_r3i1p1f1_EUROCORDEX': hists,
			   #'ECMWF-IFS-HR_r4i1p1f1_EUROCORDEX': hists,
			   #'ECMWF-IFS-HR_r5i1p1f1_EUROCORDEX': hists,
			   #'ECMWF-IFS-HR_r6i1p1f1_EUROCORDEX': hists,
			   #
			   #HighResMIP low-resolution models (LR)
			   #'EC-Earth3P_r1i1p2f1': hists,    #71km
			   #'MPI-ESM1-2-HR_r1i1p1f1': hists, #67km
			   #'CMCC-CM2-HR4_r1i1p1f1': hists,  #63km
			   #'ECMWF-IFS-LR_r1i1p1f1': hists,  #50km regridded at 100km
			   #'ECMWF-IFS-LR_r2i1p1f1': hists,
			   #'ECMWF-IFS-LR_r3i1p1f1': hists,
			   #'ECMWF-IFS-LR_r4i1p1f1': hists,
			   #'ECMWF-IFS-LR_r5i1p1f1': hists,
			   #'ECMWF-IFS-LR_r6i1p1f1': hists,
			   #'ECMWF-IFS-LR_r7i1p1f1': hists,
			   #'ECMWF-IFS-LR_r8i1p1f1': hists,
			   #'CNRM-CM6-1_r1i1p1f2': hists,      #142km
			   #'HadGEM3-GC31-LL_r1i1p1f1': hists, #135km
			   #'HadGEM3-GC31-LL_r1i2p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i3p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i4p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i5p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i6p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i7p1f1': hists,
			   #'HadGEM3-GC31-LL_r1i8p1f1': hists,
			   #
			   #Observations
			   #'obs_on_12km': hists,
			   'obs_cordex50': hists,
               'obs_cordex50_scale': hists,
                           #'GPCP': hists,
                           }
              }

# diff_btw_res = ['n512', 'n512_future',
#                 '2p2_on_n512', '2p2f_on_n512',
#                 'n512', '2p2_on_n512',
#                 'n512_future', '2p2f_on_n512', ]
# diff_btw_res=[
#             country+'_on_12km', '2p2erai_on_12km', 
#             country+'_on_12km', 'mi-ao438',
#             country+'_on_12km', '2p2eth_on_12km',
#             country+'_on_12km', '12eth_on_12km']
diff_btw_res = None

res_todo = config['histfile'].keys()  # ['ukcpobs', '2p2erai', '2p2']#config['histfile'].keys()

# ### ME: this is where you define the colour, linewidth and linestyle for each model, so add your CORDEX model names here
colordict = {#colours 'b','r','k','gold','k','y','c','m','g','aquamarine','darkgreen'
             #
	     #EUR11 models (plot in dark blue)
	     'CORDEX12': ['darkblue', 1, '-'],
	     'CCLM-MPI': ['darkblue', 1, '-'],
             'CCLM-CNRM': ['darkblue', 1, '-'],
             'CCLM-EC-EARTH': ['darkblue', 1, '-'],
             'CCLM-HadGEM': ['darkblue', 1, '-'],  
             'ALADIN-CNRM1': ['darkblue', 1, '-'],
             'ALADIN-CNRM2': ['darkblue', 1, '-'],
             'HIRHAM-CNRM': ['darkblue', 1, '-'],  
             'HIRHAM-EC-EARTH': ['darkblue', 1, '-'],			
             'HIRHAM-NorESM': ['darkblue', 1, '-'],
             'REMO-IPSL': ['darkblue', 1, '-'],
             'REMO-NorESM': ['darkblue', 1, '-'],  
             'REMO-GFDL': ['darkblue', 1, '-'],
             'WRF-IPSL': ['darkblue', 1, '-'],
             'RACMO-CNRM': ['darkblue', 1, '-'],   
             'RACMO-EC-EARTH': ['darkblue', 1, '-'],			
             'RACMO-HadGEM': ['darkblue', 1, '-'], 
             'RACMO-NorESM': ['darkblue', 1, '-'], 
             'REMO-MPI': ['darkblue', 1, '-'],
             'ALARO-CNRM': ['darkblue', 1, '-'],  
             'RCA-CNRM': ['darkblue', 1, '-'],
             'RCA-EC-EARTH': ['darkblue', 1, '-'], 
             'RCA-IPSL': ['darkblue', 1, '-'],
             'RCA-HadGEM': ['darkblue', 1, '-'],   
             'RCA-MPI': ['darkblue', 1, '-'],
             'RCA-NorESM': ['darkblue', 1, '-'],   
             #
	     #EUR44 models (plot in dark red)
             'CORDEX-44': ['darkred', 1, '-'],
             'CCLM-MPI_50km': ['darkred', 1, '-'],
             'CCLM5-CNRM_50km': ['darkred', 1, '-'],	       
             'CCLM5-EC-EARTH_50km': ['darkred', 1, '-'],	       
             'CCLM5-HadGEM_50km': ['darkred', 1, '-'],         
             'CCLM5-MIROC_50km': ['darkred', 1, '-'],	       
             'CCLM5-MPI_50km': ['darkred', 1, '-'],		       
             'ALADIN-CNRM1_50km': ['darkred', 1, '-'], 	       
             'ALADIN-CNRM2_50km': ['darkred', 1, '-'], 	       
             'HIRHAM-EC-EARTH_50km': ['darkred', 1, '-'],	       
             'WRF-IPSL_50km': ['darkred', 1, '-'],
             'WRF-CanESM_50km': ['darkred', 1, '-'],	       
             'WRF-EC-EARTH_50km': ['darkred', 1, '-'],         
             'RegCM-HadGEM_50km': ['darkred', 1, '-'], 	       
             'RACMO-EC-EARTH_50km': ['darkred', 1, '-'],	       
             'RACMO-HadGEM_50km': ['darkred', 1, '-'],         
             'REMO-MPI_50km': ['darkred', 1, '-'],
             'ALARO-CNRM_50km': ['darkred', 1, '-'],	       
             'RCA-CanESM_50km': ['darkred', 1, '-'],	       
             'RCA-CNRM_50km': ['darkred', 1, '-'],
             'RCA-CSIRO_50km': ['darkred', 1, '-'],		       
             'RCA-EC-EARTH_50km': ['darkred', 1, '-'],         
             'RCA-IPSL_50km': ['darkred', 1, '-'],
             'RCA-MIROC_50km': ['darkred', 1, '-'],		       
             'RCA-HadGEM_50km': ['darkred', 1, '-'],	       
             'RCA-MPI_50km': ['darkred', 1, '-'], 
             'RCA-NorESM_50km': ['darkred', 1, '-'],
             'RCA-GFDL_50km': ['darkred', 1, '-'],
             #
	     #CMIP5 models (plot in green)
             'CMIP5': ['chartreuse', 1, '-'],
             'MPI-ESM-LR': ['chartreuse', 1, '-'],
             'CNRM-CM5': ['chartreuse', 1, '-'],
             'EC-EARTH': ['chartreuse', 1, '-'],
             'HadGEM2-ES': ['chartreuse', 1, '-'],
             'NorESM1-M': ['chartreuse', 1, '-'],
             'IPSL-CM5A-LR': ['chartreuse', 1, '-'],
             'IPSL-CM5A-MR': ['chartreuse', 1, '-'],
             'GFDL-ESM2G': ['chartreuse', 1, '-'],
             'MIROC5': ['chartreuse', 1, '-'],
             'CanESM2': ['chartreuse', 1, '-'],
             'CSIRO-Mk3-6-0': ['chartreuse', 1, '-'],
             'GFDL-ESM2M': ['chartreuse', 1, '-'],
	     #
	     #HighResMIP HR (resolution at 50N)
	     #'PRIMAVERA': ['gold', 2, '-'],
	     'EC-Earth3P-HR_r1i1p2f1': ['gold', 2, '-'],  #36km
	     'HadGEM3-GC31-HM_r1i1p1f1': ['gold', 2, '-'],#25km
             'CMCC-CM2-VHR4_r1i1p1f1': ['gold', 2, '-'],  #18km
	     'MPI-ESM1-2-XR_r1i1p1f1': ['gold', 2, '-'],  #34km
	     'CNRM-CM6-1-HR_r1i1p1f2': ['gold', 1, '-'],  #50km
	     'ECMWF-IFS-HR_r1i1p1f1': ['gold', 2, '-'],   #25km regridded at 50km
	     'ECMWF-IFS-HR_r2i1p1f1': ['gold', 2, '-'],
	     'ECMWF-IFS-HR_r3i1p1f1': ['gold', 2, '-'],
	     'ECMWF-IFS-HR_r4i1p1f1': ['gold', 2, '-'],
	     'ECMWF-IFS-HR_r5i1p1f1': ['gold', 2, '-'],
	     'ECMWF-IFS-HR_r6i1p1f1': ['gold', 2, '-'],
	     #HighResMIP HR regridded at 50km resolution
	     'PRIMAVERA': ['orange', 1, '-'],
	     'EC-Earth3P-HR_r1i1p2f1_EUROCORDEX': ['orange', 2, '-'],
	     'HadGEM3-GC31-HM_r1i1p1f1_EUROCORDEX': ['orange', 2, '-'],
             'CMCC-CM2-VHR4_r1i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'MPI-ESM1-2-XR_r1i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'CNRM-CM6-1-HR_r1i1p1f2_EUROCORDEX': ['orange', 1, '-'],
	     'ECMWF-IFS-HR_r1i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'ECMWF-IFS-HR_r2i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'ECMWF-IFS-HR_r3i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'ECMWF-IFS-HR_r4i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'ECMWF-IFS-HR_r5i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     'ECMWF-IFS-HR_r6i1p1f1_EUROCORDEX': ['orange', 2, '-'],
	     #
	     #HighResMIP LR (resolution at 50N)
	     'PRIMAVERA LR': ['k', 1, '-'],
	     'EC-Earth3P_r1i1p2f1': ['k', 1, '-'],    #71km
	     'MPI-ESM1-2-HR_r1i1p1f1': ['k', 1, '-'], #67km	     
             'CMCC-CM2-HR4_r1i1p1f1': ['k', 1, '-'],  #64km
	     'ECMWF-IFS-LR_r1i1p1f1': ['k', 1, '-'],  #50km regridded at 100km
	     'ECMWF-IFS-LR_r2i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r3i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r4i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r5i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r6i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r7i1p1f1': ['k', 1, '-'],
	     'ECMWF-IFS-LR_r8i1p1f1': ['k', 1, '-'],
	     'CNRM-CM6-1_r1i1p1f2': ['k', 1, '-'],      #142km
	     'HadGEM3-GC31-LL_r1i1p1f1': ['k', 1, '-'], #135km
	     'HadGEM3-GC31-LL_r1i2p1f1': ['k', 1, '-'],
		'HadGEM3-GC31-LL_r1i3p1f1': ['k', 1, '-'],
	     'HadGEM3-GC31-LL_r1i4p1f1': ['k', 1, '-'],
	     'HadGEM3-GC31-LL_r1i5p1f1': ['k', 1, '-'],
	     'HadGEM3-GC31-LL_r1i6p1f1': ['k', 1, '-'],
	     'HadGEM3-GC31-LL_r1i7p1f1': ['k', 1, '-'],
	     'HadGEM3-GC31-LL_r1i8p1f1': ['k', 1, '-'],
	     #
	     #Observations
	     #'obs_on_12km': ['black', 4, '-'],
	     'obs_cordex50': ['black', 1, '-'],
         'obs_cordex50_scale': ['black', 1, '--'],
             #'GPCP': ['grey', 4, '-'],
}

# ### ME: not needed, but if you want to change the model names on the plots, enter them here.
dictnames = {'2p2f_on_n512': '2.2\,km changes',
             'n512_future': '25\,km changes',
             'n512': '25km present',
             '2p2_on_n512': '2.2\,km present', 'obs_cordex50': 'obs on CORDEX-44', 'obs_cordex50_scale': 'obs on CORDEX-44 x 1.2'}

####update11/10
sig_diff_btw_list = None# [['CORDEX-44', 'PRIMAVERA', ],
                     #['obs_cordex50', 'CORDEX-44'],
                     #['obs_cordex50', 'PRIMAVERA'], 
                     #]#'CORDEX_50']
percent_sig_per_season = {}
####endupdate11/10
for season in seasons:
    print 'season', season
    PreExt = __init__.PrecipExtremes(None, config, res_todo)

    if hist_type == '1d':
        PreExt.load_data(['masked', ], season)
        if savefig:
            # PreExt.save_subplot('masked', fractional, xlim[frequency], colordict, localdir+'/images/', season=season, bigregion=bigregion,
            #                     subregion_files=maskedreglist, nb_freq_per_day=dict_freq[frequency],
            #                     diff_btw_res=diff_btw_res,
            #                     dict_names=dictnames)
            ### ME2: change save_subplot to this, percentiles are the interval used to plot the spread based on 1000 bootstrapped sample.
            #####update11/10
            # set plot=False for more than 6 regions
            sig_val = 0.1
            percent_sig_per_season[season] = PreExt.save_subplot_with_spread('masked', fractional, xlim[frequency], colordict, localdir + '/images/', season=season,
                                bigregion=bigregion,
                                subregion_files=maskedreglist, nb_freq_per_day=dict_freq[frequency],
                                dict_names=dictnames, dict_ensemble=ensemble_dict, bootstrap_config=bootstrap_config,
                                            percentiles=[5, 95], sig_diff_btw_list=sig_diff_btw_list, #[25, 50, 75] for centiles
                                            sig_val=sig_val, vals=[[1, 10], [10, 60], [60, 400], [1, 400]], plot=True, bootstrap_or_centiles='bootstrap')

    else:
        PreExt.load_data(['masked', ], season)
        # PreExt.model_averages() #Remove if don't want average calculated from multiple models present on the graph
        if savefig:
            for reg in maskedreglist.keys():
                hist = 'hist2d_{}_{}_{:d}.json'.format(reg, wet_or_dry_spells, int(thrs))
                if wet_or_dry_spells == 'wet':
                    if frequency == 'h':
                        print reg
                        PreExt.save_2dhists(wet_or_dry_spells, thrs, 'masked', (thrs, thrs + 10), 40, season,
                                            bigregion=bigregion, subregion_files={reg: maskedreglist[reg]}, hist=hist,
                                            freq_factor=dict_freq[frequency], replace_names=False,
                                            diff_btw_res=diff_btw_res, numref_list=[], no_map =True)
                    else:
                        PreExt.save_2dhists(wet_or_dry_spells, thrs, 'masked', (
                        thrs / dict_freq[frequency], thrs / dict_freq[frequency] + 100 / dict_freq[frequency]),
                                            40 * dict_freq[frequency], season, bigregion=bigregion,
                                            subregion_files={reg: maskedreglist[reg]}, hist=hist, replace_names=True)

                else:
                    PreExt.save_2dhists(wet_or_dry_spells, thrs, 'masked', (0, 0.7 / dict_freq[frequency]),
                                        180 * dict_freq[frequency], season, bigregion=bigregion,
                                        subregion_files={reg: maskedreglist[reg]}, hist=hist,
                                        freq_factor=dict_freq[frequency], replace_names=True)

####update11/10
if 1==0:
 if percent_sig_per_season:
    print 'percent_sig_per_season', percent_sig_per_season
    subpl = __init__.plot_map(subregion_files, bigregion)
    __init__.plot_pie_inset(subregion_files, subpl, 0.8, percent_sig_per_season, localdir + '/images/', comp_col='CORDEX-44vsPRIMAVERA', comp_let=['obs_cordex50vsCORDEX-44', 'obs_cordex50vsPRIMAVERA'], plot_type=fractional, distrib=distribution, percent_interval=0.9, sig_val=int(sig_val*100))


