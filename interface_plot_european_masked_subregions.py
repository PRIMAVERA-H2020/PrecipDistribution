'''
Created on Jun 29, 2016

@author: cbates

For plotting the interface and the histogram plots to be visible
'''

import __init__
import european_masked_subregion as ems
import argparse
import numpy as np
import os
import json 

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

def join_all_files_in_one(targetfile, maskedreglist, sig_dir, end_of_file):
   dict_to_output = {'intervals':{}, 'percent_sig': {}}
   for enu, reg in enumerate(maskedreglist, 1):
      reg_file = '{}/{}_{}'.format(sig_dir, reg, end_of_file)
      with open(reg_file, 'r') as fh:
         input_dict = json.load(fh)
         #numpy_walk(input_dict
         for out in ['intervals', 'percent_sig']:
            if enu == 1:
                dict_to_output[out] = {}
            for season in ('mam', 'jja', 'son', 'djf'):
                if enu == 1:
                    dict_to_output[out][season] = {}
                dict_to_output[out][season][reg] = input_dict[out][season][reg]
   print dict_to_output
   with open(targetfile, 'w') as fh1:
      json.dump(dict_to_output, fh1, indent=2, sort_keys=True) 


def numpy_walk(node):
    assert isinstance(node, dict), 'This function only works on dictionaries'
    for key, item in node.items():
        if isinstance(item, list):
            node[key] = np.array(item)
        else:
            node[key] = item  # 0
        # numpy_walk(node)


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

#ensembles_to_compare = ['PRIMAVERA', 'CORDEX-44', 'obs_cordex50']
ensembles_to_compare = ['PRIMAVERA', 'CORDEX-11_regrid', 'obs_cordex50']
#ensembles_to_compare = [ 'PRIMAVERA']
bootstrap_or_centiles = 'centiles' 

### generate ensemble_dict: the dictionary containing all models names in each ensemble. 
### The observations are not in this as interannual variability rather than inter-member is calculated
ensemble_dict = {}
for ens in ensembles_to_compare:
    if not 'obs' in ens.lower():
        runlist, _ = ems.get_runlist_region('d', False, 'prudence', other=ens)
        ensemble_dict[ens] = runlist.keys()

## this provides the bootstrap configuration, if bootstrapping is chosen to compare the ensembles
nb_bootstrap_iter = 999999
bootstrap_config = {}
if bootstrap_or_centiles == 'bootstrap':
    bootstrap_config = {
                   'PRIMAVERA':
			{'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
			 'bootstrap_no_years': len(ensemble_dict['PRIMAVERA']),  # number of PRIMAVERA models
			 'bootstrap_sample_size': nb_bootstrap_iter,  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
 		 'bootstrap_samples_split': [0, nb_bootstrap_iter],  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
			 'bootstrap_filename': 'bootstrap_ranint_master_primavera_{}.txt'.format(nb_bootstrap_iter),  # I/O name for bootstrap years
			 'master_ref_year': 1997, # random
			},
		    #'PRIMAVERA LR':
		    #	{'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		    #	 'bootstrap_no_years': len(ensemble_dict['PRIMAVERA LR']),  # number of PRIMAVERA models
		    #	 'bootstrap_sample_size': nb_bootstrap_iter,  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    #	 'bootstrap_samples_split': [0, nb_bootstrap_iter],  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    #	 'bootstrap_filename': 'bootstrap_ranint_master_primavera_lr_{}.txt'.format(nb_bootstrap_iter),  # I/O name for bootstrap years
		    #	 'master_ref_year': 1997, # random
		    #	},
		    #'CMIP5':
		    #	{'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		    #	 'bootstrap_no_years': len(ensemble_dict['CMIP5']),  # number of PRIMAVERA models
		    #	 'bootstrap_sample_size': nb_bootstrap_iter,  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    #	 'bootstrap_samples_split': [0, nb_bootstrap_iter],  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    #	 'bootstrap_filename': 'bootstrap_ranint_master_cmip5_{}.txt'.format(nb_bootstrap_iter),  # I/O name for bootstrap years
		    #	 'master_ref_year': 1997, # random
		    #	},
		    'CORDEX-11_regrid_reduced':
		    	{'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		    	 'bootstrap_no_years': 38,#len(ensemble_dict['CORDEX12']),  # number of PRIMAVERA models
		    	 'bootstrap_sample_size': nb_bootstrap_iter,  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    	 'bootstrap_samples_split': [0, nb_bootstrap_iter],  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		    	 'bootstrap_filename': 'bootstrap_ranint_master_cordex12_{}.txt'.format(nb_bootstrap_iter),  # I/O name for bootstrap years
		    	 'master_ref_year': 1997, # random
		    	},
		    'CORDEX-44_reduced':
		       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		        'bootstrap_no_years': 26, #len(ensemble_dict['CORDEX-44_reduced']),  # number of PRIMAVERA models
		        'bootstrap_sample_size': nb_bootstrap_iter,  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		        'bootstrap_samples_split': [0, nb_bootstrap_iter],  # NB if this goes > nb_bootstrap_iter (3 digits) formatting fails
		        'bootstrap_filename': 'bootstrap_ranint_master_cordex50_{}.txt'.format(nb_bootstrap_iter),  # I/O name for bootstrap years
		        'master_ref_year': 1997, # random
		       },		       
		    #'obs_on_12km':
		    'obs_cordex50':
		       {
		       'run_IP_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
		       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
		        'bootstrap_no_years': 33,  # nb of OBS years, different by region
		        'bootstrap_sample_size': nb_bootstrap_iter,
		        'bootstrap_samples_split': [0, nb_bootstrap_iter],
		        'bootstrap_filename': 'bootstrap_ranint_master_ibp_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
		        'master_ref_year': 1971, # start year, different by region
		       },
                       'run_AL_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_alps_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },		      
                       'run_BI_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter, 
                        'bootstrap_samples_split': [0, nb_bootstrap_iter], 
                        'bootstrap_filename': 'bootstrap_ranint_master_uk_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       }, 
		       'run_CA_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_carpathian_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
		       'run_FR_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                    	'bootstrap_sample_size': nb_bootstrap_iter,
                    	'bootstrap_samples_split': [0, nb_bootstrap_iter],
                    	'bootstrap_filename': 'bootstrap_ranint_master_france_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                    	'master_ref_year': 1971, # start year, different by region
                       },
                       'run_CE_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_germany_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_MD_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_mediterranean_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_SC_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_scandinavia_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       'run_NEE_'+distribution+'.json': ###ME2 creer une entree par region qui a des observations. modifier les champs ci-dessous en fonction de ce qu'il y a dans le .json file
                       {'bootstrap_produce_file': False,  # True if want to generate bootstrapping random years rather than use exisiting file
                        'bootstrap_no_years': 35,  # nb of OBS years, different by region
                        'bootstrap_sample_size': nb_bootstrap_iter,
                        'bootstrap_samples_split': [0, nb_bootstrap_iter],
                        'bootstrap_filename': 'bootstrap_ranint_master_scandinavia_{}.txt'.format(nb_bootstrap_iter),  # change file name by region
                        'master_ref_year': 1971, # start year, different by region
                       },
                       }, 
		    }   

########## Generate bootstrapping files #############
    for ens in ensembles_to_compare:
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
			string_placeholder = "{:6d} " + multi_placeholder + "\n"
			for a in range(0, bootstrap_config[ens][region]['bootstrap_sample_size']):
			    file.write(string_placeholder.format(a, *ran_out[a]))
     else:
	    if bootstrap_config[ens]['bootstrap_produce_file']:
		with open(bootstrap_config[ens]['bootstrap_filename'], "w") as file:
		    ran_out = np.random.randint(0, bootstrap_config[ens]['bootstrap_no_years'],
						size=(bootstrap_config[ens]['bootstrap_sample_size'],
						      bootstrap_config[ens]['bootstrap_no_years']))
		    multi_placeholder = "{:2d} " * bootstrap_config[ens]['bootstrap_no_years']
		    string_placeholder = "{:6d} " + multi_placeholder + "\n"
		    for a in range(0, bootstrap_config[ens]['bootstrap_sample_size']):
			file.write(string_placeholder.format(a, *ran_out[a]))

###### x axis limits for plots #####
xlim = {'d': (0.9, 500),
        'h': (0.02, 60)}  # (1./24./4., 100)
savefig = True
localdir = '/home/users/sberthou/PrecipExtremes/'  ### must be the same as in process_data_europe, config['datadir']
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

class Europe(object):
    def __init__(self):
        self.region = Region('EUROPE',
                             lon_start=-13.0, lon_end=32.0,
                             lat_start=33.0, lat_end=71.0, #63
                             dictppt={})


_, maskedreglist = ems.get_runlist_region(frequency, n512, country, other)

country_limits = {
                  'europe': Europe().region,
                  'prudence': Europe().region,
                  'nee': Europe().region,
                  'ca': Europe().region,
                  'ibp': Europe().region,
                  'uk': Europe().region,
                  'sc': Europe().region,
                  'ce': Europe().region,
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

config = {
              'jsondir': '{}/srex_pdf_json{}'.format(localdir, json_freq[frequency]),
              'significance_dir' : '{}/sig_dir_EUR11_welch_all{}'.format(localdir, json_freq[frequency]),
              'histfile': {               #Observations
                           #'obs_on_12km': hists,
                           'obs_cordex50': hists,
                           #'obs_cordex50_scale': hists,
                          }
         }

for ens in ensemble_dict:
    for model in ensemble_dict[ens]:
        config['histfile'][model] = hists

print 'config', config

diff_btw_res = None

res_todo = config['histfile'].keys()  # ['ukcpobs', '2p2erai', '2p2']#config['histfile'].keys()

# ### ME: this is where you define the colour, linewidth and linestyle for each model, so add your CORDEX model names here
colordict = {#colours 'b','r','k','gold','k','y','c','m','g','aquamarine','darkgreen'
             #
	     #EUR11 models (plot in dark blue)
	     'CORDEX-11': ['darkblue', 1, '-'],
         'CORDEX-11_regrid': ['darkblue', 1, '-'],
	     #EUR44 models (plot in dark red)
             'CORDEX-44': ['darkred', 1, '-'],
	     #CMIP5 models (plot in green)
             'CMIP5': ['chartreuse', 1, '-'],
             'CMIP5_reduced': ['chartreuse', 1, '-'],
	     #HighResMIP HR (resolution at 50N)
	     #'PRIMAVERA': ['gold', 2, '-'],
	     #HighResMIP HR regridded at 50km resolution
	     'PRIMAVERA': ['orange', 1, '-'],
	     #HighResMIP LR (resolution at 50N)
	     'PRIMAVERA LR': ['k', 1, '-'],
	     #
	     #Observations
	     #'obs_on_12km': ['black', 4, '-'],
	     'obs_cordex50': ['black', 1, '-'],
         'obs_cordex50_scale': ['black', 1, '--'],
             #'GPCP': ['grey', 4, '-'],
}

# ### ME: not needed, but if you want to change the model names on the plots, enter them here.
dictnames = {
             'obs_cordex50': 'observations@50km', 
             'CORDEX-11_regrid': 'EUR-11@50km',
             'PRIMAVERA': 'PRIMAVERA@50km',
             'obs_cordex50_scale': 'obs on CORDEX-44_reduced x 1.2'}

### calculates the p_value of the difference for the folllowing datasets (used for the pie chart):
sig_diff_btw_list = [
                    ['obs_cordex50', 'CORDEX-11_regrid'],
                    #['obs_cordex50_scale', 'CORDEX-44'],
                    ['obs_cordex50', 'PRIMAVERA'],
                    ['CORDEX-11_regrid', 'PRIMAVERA', ]
                    #['CORDEX-44', 'PRIMAVERA', ]
                    ]
# None# [['CORDEX-44_reduced', 'PRIMAVERA', ], #should start with 'obs' in it if using terciles!
                     #['obs_cordex50', 'CORDEX-44_reduced'],
                     #['obs_cordex50', 'PRIMAVERA'], 
                     #]#'CORDEX_50']

percent_sig_per_season = {}
intervals_per_season = {}
sig_val = 0.10
plot_histograms = True
#bootstrap_or_centiles = 'centiles'  # if 'centiles': plots centiles on the histograms (25/75) and performs a t_test
                                    # if 'bootstrap': uses the bootstrap_config 
if bootstrap_or_centiles == 'centiles':
   nb_bootstrap_iter = 't_test'
   percent_to_plot = [25, 50, 75]
else:
   percent_to_plot = [5, 95]
end_of_file = '{}_{}.json'.format(int(sig_val*100), nb_bootstrap_iter)
targetfile = '{}/{}_{}'.format(config['significance_dir'], '_'.join(r for r in maskedreglist.keys()), end_of_file)

if not plot_histograms:
  if not os.path.exists(targetfile):
    if os.path.exists('{}/{}_{}'.format(config['significance_dir'], maskedreglist.keys()[0], end_of_file)):
        print 'joining all regional files together into {}'.format(targetfile)
        join_all_files_in_one(targetfile, maskedreglist.keys(), config['significance_dir'], end_of_file)

if not os.path.exists(targetfile):
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
            ### change save_subplot to save_subplot_with_spread, 
            # percentiles are the interval used to plot the spread based on 1000 bootstrapped sample.
            # set plot=False for more than 6 regions
            percent_sig_per_season[season], intervals_per_season[season] = PreExt.save_subplot_with_spread('masked', 
                                fractional, # plot type
                                xlim[frequency], #x limits for the graph
                                colordict, #defines the colour of the line
                                localdir + '/images/', 
                                season=season,
                                bigregion=bigregion,
                                subregion_files=maskedreglist, 
                                nb_freq_per_day=dict_freq[frequency],
                                dict_names=dictnames, 
                                dict_ensemble=ensemble_dict, 
                                bootstrap_config=bootstrap_config,
                                percentiles=percent_to_plot, # [25, 50, 75] for centiles
                                sig_diff_btw_list=sig_diff_btw_list, # will be used as input for the pie chart
                                sig_val=sig_val, # threshold of p_value to define "significantly different" 
                                vals='terciles', #[[1, 10], [10, 60], [60, 400], [1, 400]], # precipitation intervals on which to calculate significant difference
                                plot=plot_histograms, # set to False if you only want to calculate significance for the pie plot
                                bootstrap_or_centiles=bootstrap_or_centiles) # can use either centiles over the ensemble or bootstrap spread to define the significance envelope.


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

 if percent_sig_per_season:
      if not os.path.isdir(config['significance_dir']):
          os.makedirs(config['significance_dir'])
      with open(targetfile, 'w') as fh:
          output = {'percent_sig': percent_sig_per_season, 'intervals': intervals_per_season}
          json.dump(output, fh, indent=2, sort_keys=True)


####update11/10
# this is where the pie plot is drawn:
else:
 print('{} exists, drawing the pie charts'.format(targetfile))
 with open(targetfile, 'r') as fh:
   input_dict = json.load(fh)
   percent_sig_per_season = input_dict['percent_sig']
   if percent_sig_per_season:
    
    print 'percent_sig_per_season', percent_sig_per_season
    ### 1) plot the map with masked subregions
    subpl = __init__.plot_map(subregion_files, bigregion)
    ### 2) plot the pie inset for each subregion
    __init__.plot_pie_inset(subregion_files, 
                            subpl, 
                            0.8, # size of the inset pie - not important 
                            percent_sig_per_season, # this was calculated by save_subplot_with_spread above
                            localdir + '/images/', 
                            comp_col='CORDEX-11_regridvsPRIMAVERA', # which datasets to compare (needs to be one defined in sig_diff_btw_list)
                            comp_let = ['obs_cordex50vsCORDEX-11_regrid', 'obs_cordex50vsPRIMAVERA'], #['obs_cordex50vsCORDEX-11_regrid_reduced', 'obs_cordex50vsPRIMAVERA'], # score of each dataset against obs to be compared, defined in sig_diff_btw_list
                            plot_type=fractional, # contrib
                            distrib=distribution, # exponential100 
                            percent_interval=0.9, # on how much of the interval are the datasets significantly different
                            sig_val=int(sig_val*100)) # what is the value meaning "significantly different
