'''
Created on Sep 13, 2016

@author: sberthou
'''
import subprocess
import collections

projdir = '/project/hires_rcm/hadek/'
datadir = '/home/users/sberthou/PrecipExtremes/'  ### change path to where your PRUDENCE mask is
primavera_path = '/gws/nopw/j04/primavera5/stream1/CMIP6/HighResMIP/{}/{}/hist-1950/{}/day/pr/g?/{}/pr_day_{}_hist-1950_{}_g?_*.nc'


def get_runlist_region(frequency, n512, country, other=''):
    '''
   This is to centralised the definition of masks and simulations and to allow to choose what runlist and which mask to use according to four parameters:

   :param str frequency: 'h', 'd'
   :param bool n512: to use the regridded versions on n512
   :param str other: if 'ukcpm_runs': to use Giorgia's runs; or 'ammacatch' 
   :param str country: string of the country to subset
   :return: runlist and masks of the region
   '''

    if n512:
        if other == 'ammacatch':
            if frequency == '15mn':
                cty_or_n512 = 'u-ad251'
            if frequency == '3h':
                cty_or_n512 = 'TRMM'
        else:
            cty_or_n512 = 'n512'

    else:
        if other == 'ammacatch':
            cty_or_n512 = 'AMMA-CATCH'
        else:
            cty_or_n512 = 'mi-ao438'  # country

    dict_subregions = {'prudence': collections.OrderedDict([
                                   # ('NEE', datadir+'mask_PRUDENCEreg_NEE.nc'),
                                    ('IP', datadir+'mask_PRUDENCEreg_IP.nc'), 
                                    ('FR', datadir+'mask_PRUDENCEreg_FR.nc'),
                                   # ('MD', datadir+'mask_PRUDENCEreg_MD.nc'),
                                    ('CA', datadir+'mask_PRUDENCEreg_CA.nc'),
                                   # ('CE', datadir+'mask_PRUDENCEreg_CE.nc'),
                                   # ('SC', datadir+'mask_PRUDENCEreg_SC.nc'),
                                    ('AL', datadir+'mask_PRUDENCEreg_AL.nc'),
                                    ##'EE': datadir+'mask_PRUDENCEreg_EE.nc',
                                    ('BI', datadir+'mask_PRUDENCEreg_BI.nc'),
                                    ###'CH': datadir+'mask_PRUDENCEreg_CH.nc',
                                    ### 'TCZ': datadir+'mask_PRUDENCEreg_TCZ.nc'
                                    ]),
                        'uk': {'BI': datadir+'mask_PRUDENCEreg_BI.nc',},
                        'ca': {'CA': datadir+'mask_PRUDENCEreg_CA.nc',},
                        'alps': {'AL': datadir+'mask_PRUDENCEreg_AL.nc',},
                        'ibp': {'IP': datadir+'mask_PRUDENCEreg_IP.nc',},
                        'france': {'FR': datadir+'mask_PRUDENCEreg_FR.nc',},
                        'nee': {'NEE': datadir+'mask_PRUDENCEreg_NEE.nc',},
                        'md': {'MD': datadir+'mask_PRUDENCEreg_MD.nc',},
                        'ce': {'CE': datadir+'mask_PRUDENCEreg_CE.nc'},
                        'sc': {'SC': datadir+'mask_PRUDENCEreg_SC.nc'},

                       }

    maskedreglist = dict_subregions[country.lower()]

    runlist_frequency = {'d':
                             {
                              'alps_cordex50': '/gws/nopw/j04/primavera3/cache/sberthou/alps_regrid_CORDEX50/RapdD_al05_CORDEX50_*.nc',
                              'ibp_cordex50': '/gws/nopw/j04/primavera3/cache/sberthou/IBP_regrid_CORDEX50/*.nc',
                              'uk_cordex50': '/gws/nopw/j04/primavera3/cache/sberthou/ukcpobs_regrid_CORDEX50/*.nc',
                              'ca_cordex50': '/gws/nopw/j04/primavera3/cache/sberthou/CARPATCLIM_regrid_CORDEX50/CCLIM_pre19612010_EUROCORDEX.nc',
                              'france_cordex50': '/gws/nopw/j04/primavera3/cache/sberthou/FRANCE_regrid_CORDEX50/ForcPRCPT*.nc',
                              },

                         'h':
                             {
                              },
                         '15mn':
                             {},
                         '3h': {
                         }
                         }

    if n512:
        if other == 'ammacatch':
            runlist = {'cp4': runlist_frequency[frequency]['cp4']}
        else:
            runlist = {  # country+'_on_n512':runlist_frequency[frequency][country+'_on_n512'],
                #    '2p2erai_on_n512':runlist_frequency[frequency]['2p2erai_on_n512'],
                '2p2_on_n512': runlist_frequency[frequency]['2p2_on_n512'],
                'n512': runlist_frequency[frequency]['n512'],
                '2p2f_on_n512': runlist_frequency[frequency]['2p2f_on_n512'],
                'n512_future': runlist_frequency[frequency]['n512_future'],
            }

    elif other == 'CORDEX-44':
        cordex_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50/'
        runlist = {
             #EUR44 models
                     'CCLM-MPI_50km': cordex_dir+'/pr/CLMcom/MPI-M-MPI-ESM-LR/CLMcom-CCLM4-8-17/r1i1p1/*_v1_*.nc',
                     'CCLM5-CNRM_50km': cordex_dir+'/pr/ETH/CNRM-CERFACS-CNRM-CM5/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     'CCLM5-EC-EARTH_50km': cordex_dir+'/pr/ETH/ICHEC-EC-EARTH/CLMcom-CCLM5-0-6/r12i1p1/*_v1_*.nc',
                     'CCLM5-HadGEM_50km': cordex_dir+'/pr/ETH/MOHC-HadGEM2-ES/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc', #problem? -> check
                     'CCLM5-MIROC_50km': cordex_dir+'/pr/ETH/MIROC-MIROC5/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     'CCLM5-MPI_50km': cordex_dir+'/pr/ETH/MPI-M-MPI-ESM-LR/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     'HIRHAM-EC-EARTH_50km': cordex_dir+'/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v1_*.nc',
                     'WRF-IPSL_50km': cordex_dir+'/pr/IPSL-INERIS/IPSL-IPSL-CM5A-MR/IPSL-INERIS-WRF331F/r1i1p1/*_v1_*.nc', #problems
                     'WRF-EC-EARTH_50km': cordex_dir+'/pr/NUIM/ICHEC-EC-EARTH/NUIM-WRF341E/r1i1p1/*_v1_*.nc', #problems!!!!check it worked
                     'RACMO-EC-EARTH_50km': cordex_dir+'/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r1i1p1/*_v1_*.nc',
                     #'RACMO-EC-EARTH12_50km': cordex_dir+'/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r12i1p1/*_v1_*.nc', ######
                     'RACMO-HadGEM_50km': cordex_dir+'/pr/KNMI/MOHC-HadGEM2-ES/KNMI-RACMO22E/r1i1p1/*_v2_*.nc',
                     'REMO-MPI_50km': cordex_dir+'/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r1i1p1/*_v1_*.nc',
                     #'REMO-MPI2_50km': cordex_dir+'/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r2i1p1/*_v1_*.nc',#####
                     'RCA-CanESM_50km': cordex_dir+'/pr/SMHI/CCCma-CanESM2/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-CNRM_50km': cordex_dir+'/pr/SMHI/CNRM-CERFACS-CNRM-CM5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-CSIRO_50km': cordex_dir+'/pr/SMHI/CSIRO-QCCCE-CSIRO-Mk3-6-0/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-EC-EARTH_50km': cordex_dir+'/pr/SMHI/ICHEC-EC-EARTH/SMHI-RCA4/r12i1p1/*_v1_*.nc',
                     'RCA-IPSL_50km': cordex_dir+'/pr/SMHI/IPSL-IPSL-CM5A-MR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-MIROC_50km': cordex_dir+'/pr/SMHI/MIROC-MIROC5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-HadGEM_50km': cordex_dir+'/pr/SMHI/MOHC-HadGEM2-ES/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-MPI_50km': cordex_dir+'/pr/SMHI/MPI-M-MPI-ESM-LR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-NorESM_50km': cordex_dir+'/pr/SMHI/NCC-NorESM1-M/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-GFDL_50km': cordex_dir+'/pr/SMHI/NOAA-GFDL-GFDL-ESM2M/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'ALADIN-CNRM1_50km': cordex_dir+'/pr/CNRM/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN53/r1i1p1/*_v1_*.nc',
                     #'ALADIN-CNRM2_50km': cordex_dir+'/pr/HMS/CNRM-CERFACS-CNRM-CM5/HMS-ALADIN52/r1i1p1/*_v1_*.nc',
                     'WRF-CanESM_50km': cordex_dir+'/pr/UCAN/CCCma-CanESM2/UCAN-WRF341I/r1i1p1/*_v2_*.nc', ####
                     'RegCM-HadGEM_50km': cordex_dir+'/pr/ICTP/MOHC-HadGEM2-ES/ICTP-RegCM4-3/r1i1p1/*_v1_*.nc', ####
                     'ALARO-CNRM_50km': cordex_dir+'/pr/RMIB-UGent/CNRM-CERFACS-CNRM-CM5/RMIB-UGent-ALARO-0/r1i1p1/*_v1_*.nc',
                   }
    elif other == 'CORDEX-44_reduced':
        cordex_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50/'
        runlist = {
             #EUR44 models
                     'CCLM-MPI_50km': cordex_dir+'/pr/CLMcom/MPI-M-MPI-ESM-LR/CLMcom-CCLM4-8-17/r1i1p1/*_v1_*.nc',
                     'CCLM5-CNRM_50km': cordex_dir+'/pr/ETH/CNRM-CERFACS-CNRM-CM5/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     'CCLM5-EC-EARTH_50km': cordex_dir+'/pr/ETH/ICHEC-EC-EARTH/CLMcom-CCLM5-0-6/r12i1p1/*_v1_*.nc',
                     'CCLM5-HadGEM_50km': cordex_dir+'/pr/ETH/MOHC-HadGEM2-ES/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc', #problem? -> check
                     'CCLM5-MPI_50km': cordex_dir+'/pr/ETH/MPI-M-MPI-ESM-LR/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     'HIRHAM-EC-EARTH_50km': cordex_dir+'/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v1_*.nc',
                     'WRF-EC-EARTH_50km': cordex_dir+'/pr/NUIM/ICHEC-EC-EARTH/NUIM-WRF341E/r1i1p1/*_v1_*.nc', #problems!!!!check it worked
                     'RACMO-EC-EARTH_50km': cordex_dir+'/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r1i1p1/*_v1_*.nc',
                     'RACMO-HadGEM_50km': cordex_dir+'/pr/KNMI/MOHC-HadGEM2-ES/KNMI-RACMO22E/r1i1p1/*_v2_*.nc',
                     'REMO-MPI_50km': cordex_dir+'/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r1i1p1/*_v1_*.nc',
                     'RCA-CNRM_50km': cordex_dir+'/pr/SMHI/CNRM-CERFACS-CNRM-CM5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-EC-EARTH_50km': cordex_dir+'/pr/SMHI/ICHEC-EC-EARTH/SMHI-RCA4/r12i1p1/*_v1_*.nc',
                     'RCA-HadGEM_50km': cordex_dir+'/pr/SMHI/MOHC-HadGEM2-ES/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-MPI_50km': cordex_dir+'/pr/SMHI/MPI-M-MPI-ESM-LR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'ALADIN-CNRM1_50km': cordex_dir+'/pr/CNRM/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN53/r1i1p1/*_v1_*.nc',
                     'RegCM-HadGEM_50km': cordex_dir+'/pr/ICTP/MOHC-HadGEM2-ES/ICTP-RegCM4-3/r1i1p1/*_v1_*.nc', ####
                     'ALARO-CNRM_50km': cordex_dir+'/pr/RMIB-UGent/CNRM-CERFACS-CNRM-CM5/RMIB-UGent-ALARO-0/r1i1p1/*_v1_*.nc',
                   }

    elif other == 'CORDEX-11':
        cordex11_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/'
        runlist = {
             #EUR11 models on native grid
                     'CCLM-MPI': cordex11_dir+'/EUR-11/historical/day/pr/CLMcom/MPI-M-MPI-ESM-LR/CLMcom-CCLM4-8-17/r1i1p1/*.nc',
                     'CCLM-CNRM': cordex11_dir+'/EUR-11/historical/day/pr/CLMcom/CNRM-CERFACS-CNRM-CM5/CLMcom-CCLM4-8-17/r1i1p1/*.nc',
                     #'CCLM-EC-EARTH': cordex11_dir+'/EUR-11/historical/day/pr/CLMcom/ICHEC-EC-EARTH/CLMcom-CCLM4-8-17/r12i1p1/*.nc',
                     'CCLM-HadGEM': cordex11_dir+'/EUR-11/historical/day/pr/CLMcom/MOHC-HadGEM2-ES/CLMcom-CCLM4-8-17/r1i1p1/*.nc',

                     'CCLM-ETH-MPI': cordex11_dir+'EUR-11/historical/day/pr/CLMcom-ETH/MPI-M-MPI-ESM-LR/CLMcom-ETH-COSMO-crCLIM-v1-1/r1i1p1/*.nc', ## new, r2 and r3 also exist
                     'CCLM-ETH-EC-Earth': cordex11_dir+'EUR-11/historical/day/pr/CLMcom-ETH/ICHEC-EC-EARTH/CLMcom-ETH-COSMO-crCLIM-v1-1/r12i1p1/*.nc',## new
                     'CCLM-ETH-NorESM': cordex11_dir+'EUR-11/historical/day/pr/CLMcom-ETH/NCC-NorESM1-M/CLMcom-ETH-COSMO-crCLIM-v1-1/r1i1p1/*.nc',## new

                     'HIRHAM-CNRM1': cordex11_dir+'/EUR-11/historical/day/pr/DMI/CNRM-CERFACS-CNRM-CM5/DMI-HIRHAM5/r1i1p1/*_v1_*.nc',
                     #'HIRHAM-CNRM2': cordex11_dir+'/EUR-11/historical/day/pr/DMI/CNRM-CERFACS-CNRM-CM5/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',
                     'HIRHAM-EC-EARTH1': cordex11_dir+'/EUR-11/historical/day/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r1i1p1/*_v1_*.nc',
                     #'HIRHAM-EC-EARTH12': cordex11_dir+'/EUR-11/historical/day/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r12i1p1/*_v1_*.nc',
                     #'HIRHAM-EC-EARTH31': cordex11_dir+'/EUR-11/historical/day/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v1_*.nc',
                     #'HIRHAM-EC-EARTH32': cordex11_dir+'/EUR-11/historical/day/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v2_*.nc',
                     'HIRHAM-HadGEM1': cordex11_dir+'/EUR-11/historical/day/pr/DMI/MOHC-HadGEM2-ES/DMI-HIRHAM5/r1i1p1/*_v1_*.nc',
                     #'HIRHAM-HadGEM2': cordex11_dir+'/EUR-11/historical/day/pr/DMI/MOHC-HadGEM2-ES/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',
                     'HIRHAM-NorESM2': cordex11_dir+'/EUR-11/historical/day/pr/DMI/NCC-NorESM1-M/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',
                     'HIRAM-MPI': cordex11_dir+'EUR-11/historical/day/pr/DMI/MPI-M-MPI-ESM-LR/DMI-HIRHAM5/r1i1p1/*.nc',              ## new
                     #'HIRHAM-NorESM3': cordex11_dir+'/EUR-11/historical/day/pr/DMI/NCC-NorESM1-M/DMI-HIRHAM5/r1i1p1/*_v3_*.nc',

                     'REMO-IPSL': cordex11_dir+'/EUR-11/historical/day/pr/GERICS/IPSL-IPSL-CM5A-LR/GERICS-REMO2015/r1i1p1/*_v1_*.nc',
                     'REMO-NorESM': cordex11_dir+'/EUR-11/historical/day/pr/GERICS/NCC-NorESM1-M/GERICS-REMO2015/r1i1p1/*_v1_*.nc',
                     'REMO-GFDL': cordex11_dir+'/EUR-11/historical/day/pr/GERICS/NOAA-GFDL-GFDL-ESM2G/GERICS-REMO2015/r1i1p1/*_v1_*.nc',
                     'REMO-MPI': cordex11_dir+'/EUR-11/historical/day/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r1i1p1/*_v1_*.nc',
                     #'REMO-MPI2': cordex11_dir+'/EUR-11/historical/day/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r2i1p1/*_v1_*.nc',

                     ##'WRF-IPSL': cordex11_dir+'/EUR-11/historical/day/pr/IPSL-INERIS/IPSL-IPSL-CM5A-MR/IPSL-INERIS-WRF331F/r1i1p1/*_v1_*.nc',

                     'RACMO-CNRM': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/CNRM-CERFACS-CNRM-CM5/KNMI-RACMO22E/r1i1p1/*_v2_*.nc',
                     ##'RACMO-EC-EARTH': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r1i1p1/*_v1_*.nc',
                     'RACMO-EC-EARTH12': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r12i1p1/*_v1_*.nc',
                     #'RACMO-EC-EARTH3': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r3i1p1/*_v1_*.nc',
                     'RACMO-HadGEM': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/MOHC-HadGEM2-ES/KNMI-RACMO22E/r1i1p1/*_v2_*.nc',
                     'RACMO-NorESM': cordex11_dir+'/EUR-11/historical/day/pr/KNMI/NCC-NorESM1-M/KNMI-RACMO22E/r1i1p1/*_v1_*.nc',
                     'RACMO-IPSL': cordex11_dir + '/EUR-11/historical/day/pr/KNMI/IPSL-IPSL-CM5A-MR/KNMI-RACMO22E/r1i1p1/*.nc', ## new
                     'RACMO-MPI': cordex11_dir + '/EUR-11/historical/day/pr/KNMI/MPI-M-MPI-ESM-LR/KNMI-RACMO22E/r1i1p1/*.nc', ## new

                     'RCA-CNRM': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/CNRM-CERFACS-CNRM-CM5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-EC-EARTH': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/ICHEC-EC-EARTH/SMHI-RCA4/r12i1p1/*_v1_*.nc',
                     #'RCA-EC-EARTH3': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/ICHEC-EC-EARTH/SMHI-RCA4/r3i1p1/*_v1_*.nc',
                     'RCA-IPSL': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/IPSL-IPSL-CM5A-MR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-HadGEM': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/MOHC-HadGEM2-ES/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     'RCA-MPI': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/MPI-M-MPI-ESM-LR/SMHI-RCA4/r1i1p1/*_v1a_*.nc',
                     'RCA-NorESM': cordex11_dir+'/EUR-11/historical/day/pr/SMHI/NCC-NorESM1-M/SMHI-RCA4/r1i1p1/*_v1_*.nc',

                     ##'ALADIN-CNRM1': cordex11_dir+'/EUR-11/historical/day/pr/CNRM/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN53/r1i1p1/*_v1_*.nc',
                     ##'ALADIN-CNRM2': cordex11_dir+'/EUR-11/historical/day/pr/CNRM/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/r1i1p1/*_v2_*.nc',
                     ##'ALARO-CNRM': cordex11_dir+'/EUR-11/historical/day/pr/RMIB-UGent/CNRM-CERFACS-CNRM-CM5/RMIB-UGent-ALARO-0/r1i1p1/*_v1_*.nc',
                     #'ALADIN-HadGEM2': cordex11_dir+'EUR-11/historical/day/pr/CNRM/MOHC-HadGEM2-ES/CNRM-ALADIN63/r1i1p1/*.nc', ## new
                     #'ALADIN-MPI': cordex11_dir+'EUR-11/historical/day/pr/CNRM/MPI-M-MPI-ESM-LR/CNRM-ALADIN63/r1i1p1/*.nc', ##new

                     'MOHC_ICHEC-EC-EARTH': cordex11_dir+'EUR-11/historical/day/pr/MOHC/ICHEC-EC-EARTH/MOHC-HadREM3-GA7-05/r12i1p1/*.nc', ##new
                     'MOHC_HadGEM2': cordex11_dir+'EUR-11/historical/day/pr/MOHC/MOHC-HadGEM2-ES/MOHC-HadREM3-GA7-05/r1i1p1/*.nc', ##new

                     'WRF-CNRM': cordex11_dir+'EUR-11/historical/day/pr/IPSL/CNRM-CERFACS-CNRM-CM5/IPSL-WRF381P/r1i1p1/*.nc',  ##new 
                     #'WRF-IPSL': cordex11_dir+'EUR-11/historical/day/pr/IPSL/IPSL-IPSL-CM5A-MR/IPSL-WRF381P/r1i1p1/*.nc',      ##new
                     'WRF-NorESM': cordex11_dir+ 'EUR-11/historical/day/pr/IPSL/NCC-NorESM1-M/IPSL-WRF381P/r1i1p1/*.nc',     #new
                     'WRF-HadGEM2': cordex11_dir+ 'EUR-11/historical/day/pr/IPSL/MOHC-HadGEM2-ES/IPSL-WRF381P/r1i1p1/*.nc',    #new
           }
    elif other == 'CORDEX-11_regrid':
            cordex11_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/regridding_cordex12/'
            runlist = {
             'CCLM-MPI_regrid': cordex11_dir+'/CLMcom/*_MPI-M-MPI-ESM-LR*r1i1p1*.nc',
             'CCLM-CNRM_regrid': cordex11_dir+'/CLMcom/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*.nc',
             'CCLM-EC-EARTH_regrid': cordex11_dir+'/CLMcom/*_ICHEC-EC-EARTH*r12i1p1*.nc',
             'CCLM-HadGEM_regrid': cordex11_dir+'/CLMcom/*_MOHC-HadGEM2-ES*r1i1p1*.nc',

             'CCLM-ETH-MPI_regrid': cordex11_dir+'/CLMcom-ETH/pr_EUR-11_MPI-M-MPI-ESM-LR_historical*.nc', ## new, r2 and r3 also exist
             'CCLM-ETH-EC-Earth_regrid': cordex11_dir+'/CLMcom-ETH/pr_EUR-11_ICHEC-EC-EARTH_historical*.nc',## new
             'CCLM-ETH-NorESM_regrid': cordex11_dir+'/CLMcom-ETH/pr_EUR-11_NCC-NorESM1-M_historical*.nc',## new

             'HIRHAM-CNRM1_regrid': cordex11_dir+'/DMI/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*_v1_*.nc',
             'HIRHAM-CNRM2_regrid': cordex11_dir+'/DMI/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*_v2_*.nc',
             'HIRHAM-EC-EARTH1_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH*r1i1p1*_v1_*.nc',
             ##'HIRHAM-EC-EARTH12_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH/DMI-HIRHAM5/r12i1p1/*_v1_*.nc',
             ##'HIRHAM-EC-EARTH31_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v1_*.nc',
             'HIRHAM-EC-EARTH32_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v2_*.nc',
             'HIRHAM-HadGEM1_regrid': cordex11_dir+'/DMI/*_MOHC-HadGEM2-ES*_v1_*.nc',
             'HIRHAM-HadGEM2_regrid': cordex11_dir+'/DMI/*_MOHC-HadGEM2-ES/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',
             'HIRHAM-NorESM2_regrid': cordex11_dir+'/DMI/*_NCC-NorESM1-M/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',
             'HIRHAM-NorESM3_regrid': cordex11_dir+'/DMI/*_NCC-NorESM1-M*_v3_*.nc',
             'HIRAM-MPI_regrid': cordex11_dir + 'DMI_2/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_*.nc',

             'REMO-IPSL_regrid': cordex11_dir+'/GERICS/*_IPSL-IPSL-CM5A-LR*_v1_*.nc',
             'REMO-NorESM_regrid': cordex11_dir+'/GERICS/*_NCC-NorESM1-M*_v1_*.nc',
             'REMO-GFDL_regrid': cordex11_dir+'/GERICS/*_NOAA-GFDL-GFDL-ESM2G*_v1_*.nc',
             'REMO-MPI2_regrid': cordex11_dir+'/GERICS_2/*_MPI-M-MPI-ESM-LR*r3i1p1*v1_*.nc',
                    ##'REMO-MPI_regrid': cordex11_dir+'/MPI-CSC/*_MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r1i1p1/*_v1_*.nc',

                     'RACMO-CNRM_regrid': cordex11_dir+'/KNMI/*_CNRM-CERFACS-CNRM-CM*_v2_*.nc',
                     'RACMO-EC-EARTH_regrid': cordex11_dir+'/KNMI/*_ICHEC-EC-EARTH*r1i1p1*_v1_*.nc',
             ##'RACMO-EC-EARTH12_regrid': cordex11_dir+'/KNMI/*_ICHEC-EC-EARTH*_v1_*.nc',
             ##'RACMO-EC-EARTH3_regrid': cordex11_dir+'/KNMI/*_ICHEC-EC-EARTH/KNMI-RACMO22E/r3i1p1/*_v1_*.nc',
                     'RACMO-HadGEM_regrid': cordex11_dir+'/KNMI/*_MOHC-HadGEM2-ES*_v2_*.nc',
                     'RACMO-NorESM_regrid': cordex11_dir+'/KNMI/*_NCC-NorESM1-M*_v1_*.nc',
                     'RACMO-IPSL_regrid': cordex11_dir + '/KNMI_2/*IPSL-CM5A-MR*r1i1p1*.nc', ## new
                     'RACMO-MPI_regrid': cordex11_dir + '/KNMI_2/*MPI-M-MPI-ESM-LR*r1i1p1*.nc', ## new

                     'RCA-CNRM_regrid': cordex11_dir+'/SMHI/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'RCA-EC-EARTH_regrid': cordex11_dir+'/SMHI/*_ICHEC-EC-EARTH*r3i1p1*_v1_*.nc',
             ##'RCA-EC-EARTH3_regrid': cordex11_dir+'/SMHI/*_ICHEC-EC-EARTH/SMHI-RCA4/r3i1p1/*_v1_*.nc',
                     'RCA-IPSL_regrid': cordex11_dir+'/SMHI/*_IPSL-IPSL-CM5A-MR*_v1_*.nc',
                     'RCA-HadGEM_regrid': cordex11_dir+'/SMHI/*_MOHC-HadGEM2-ES*_v1_*.nc',
                     'RCA-MPI_regrid': cordex11_dir+'/SMHI/*_MPI-M-MPI-ESM-LR*_v1a_*.nc',
                     'RCA-NorESM_regrid': cordex11_dir+'/SMHI/*_NCC-NorESM1-M*_v1_*.nc',

                     'ALADIN-CNRM1_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'ALADIN-CNRM2_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/r1i1p1/*_v2_*.nc',
                     'ALARO-CNRM_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'ALADIN-HadGEM2_regrid': cordex11_dir+'CNRM_2/pr_EUR-11_MOHC-HadGEM2-ES_historical_*.nc', ##new
                     'ALADIN-MPI_regrid': cordex11_dir + 'CNRM_2/pr_EUR-11_MPI-M-MPI-ESM-LR_historical*.nc', ##new
                     
                     'WRF-CNRM_regrid': cordex11_dir+'IPSL/*CNRM-CERFACS-CNRM-CM5*r1i1p1*.nc',  ##new 
                     'WRF-IPSL_regrid': cordex11_dir+'IPSL/*IPSL-IPSL-CM5A-MR*r1i1p1*.nc',      ##new
                     'WRF-NorESM_regrid': cordex11_dir+ 'IPSL/*NCC-NorESM1-M*r1i1p1*.nc',     #new
                     'WRF-HadGEM2_regrid': cordex11_dir+ 'IPSL/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_IPSL-WRF381P_v1*.nc', # new

                     'MOHC_EC-EARTH_regrid': cordex11_dir+'MOHC/*ICHEC-EC-EARTH*r12i1p1*.nc', ##new
                     'MOHC_HadGEM2_regrid': cordex11_dir+'MOHC/*MOHC-HadGEM2-ES*r1i1p1*.nc', ##new
                     
                     'RegCM-MPI_regrid': cordex11_dir + 'ICTP/*MPI*.nc', ## new
                     'RegCM_HadGEM2_regrid': cordex11_dir + 'ICTP/*HadGEM*.nc', ##new
                     # 'RegCM-MPI': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP/MPI-M-MPI-ESM-LR/ICTP-RegCM4-6/r1i1p1/pr*.nc',
                     # 'RegCM-HadGEM2': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP//MOHC-HadGEM2-ES/ICTP-RegCM4-6/r1i1p1/p*',

            }
    elif other == 'CORDEX-11_regrid_reduced':
            cordex11_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/regridding_cordex12/'
            runlist = {
             'CCLM-MPI_regrid': cordex11_dir+'/CLMcom/*_MPI-M-MPI-ESM-LR*r1i1p1*.nc',
             'CCLM-CNRM_regrid': cordex11_dir+'/CLMcom/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*.nc',
             'CCLM-EC-EARTH_regrid': cordex11_dir+'/CLMcom/*_ICHEC-EC-EARTH*r12i1p1*.nc',
             'CCLM-HadGEM_regrid': cordex11_dir+'/CLMcom/*_MOHC-HadGEM2-ES*r1i1p1*.nc',

             'CCLM-ETH-MPI_regrid': cordex11_dir+'/CLMcom-ETH/pr_EUR-11_MPI-M-MPI-ESM-LR_historical*.nc', ## new, r2 and r3 also exist
             'CCLM-ETH-EC-Earth_regrid': cordex11_dir+'/CLMcom-ETH/pr_EUR-11_ICHEC-EC-EARTH_historical*.nc',## new

             'HIRHAM-CNRM1_regrid': cordex11_dir+'/DMI/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*_v1_*.nc',
             'HIRHAM-EC-EARTH1_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH*r1i1p1*_v1_*.nc',
             'HIRHAM-HadGEM1_regrid': cordex11_dir+'/DMI/*_MOHC-HadGEM2-ES*_v1_*.nc',
             'HIRAM-MPI_regrid': cordex11_dir + 'DMI_2/pr_EUR-11_MPI-M-MPI-ESM-LR_historical_*.nc',
             'HIRHAM-CNRM2_regrid': cordex11_dir+'/DMI/*_CNRM-CERFACS-CNRM-CM5*r1i1p1*_v2_*.nc',
             'HIRHAM-EC-EARTH32_regrid': cordex11_dir+'/DMI/*_ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v2_*.nc',
             'HIRHAM-HadGEM2_regrid': cordex11_dir+'/DMI/*_MOHC-HadGEM2-ES/DMI-HIRHAM5/r1i1p1/*_v2_*.nc',

             'REMO-MPI2_regrid': cordex11_dir+'/GERICS_2/*_MPI-M-MPI-ESM-LR*r3i1p1*v1_*.nc',

                     'RACMO-CNRM_regrid': cordex11_dir+'/KNMI/*_CNRM-CERFACS-CNRM-CM*_v2_*.nc',
                     'RACMO-EC-EARTH_regrid': cordex11_dir+'/KNMI/*_ICHEC-EC-EARTH*r1i1p1*_v1_*.nc',
                     'RACMO-HadGEM_regrid': cordex11_dir+'/KNMI/*_MOHC-HadGEM2-ES*_v2_*.nc',
                     'RACMO-MPI_regrid': cordex11_dir + '/KNMI_2/*MPI-M-MPI-ESM-LR*r1i1p1*.nc', ## new

                     'RCA-CNRM_regrid': cordex11_dir+'/SMHI/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'RCA-EC-EARTH_regrid': cordex11_dir+'/SMHI/*_ICHEC-EC-EARTH*r3i1p1*_v1_*.nc',
                     'RCA-HadGEM_regrid': cordex11_dir+'/SMHI/*_MOHC-HadGEM2-ES*_v1_*.nc',
                     'RCA-MPI_regrid': cordex11_dir+'/SMHI/*_MPI-M-MPI-ESM-LR*_v1a_*.nc',

                     'ALADIN-CNRM1_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'ALARO-CNRM_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5*_v1_*.nc',
                     'ALADIN-CNRM2_regrid': cordex11_dir+'/CNRM/*_CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/r1i1p1/*_v2_*.nc', 
                     'ALADIN-HadGEM2_regrid': cordex11_dir+'CNRM_2/pr_EUR-11_MOHC-HadGEM2-ES_historical_*.nc', ##new
                     'ALADIN-MPI_regrid': cordex11_dir + 'CNRM_2/pr_EUR-11_MPI-M-MPI-ESM-LR_historical*.nc', ##new

                     'WRF-CNRM_regrid': cordex11_dir+'IPSL/*CNRM-CERFACS-CNRM-CM5*r1i1p1*.nc',  ##new 
                     'WRF-HadGEM2_regrid': cordex11_dir+ 'IPSL/pr_EUR-11_MOHC-HadGEM2-ES_historical_r1i1p1_IPSL-WRF381P_v1*.nc', # new

                     'MOHC_EC-EARTH_regrid': cordex11_dir+'MOHC/*ICHEC-EC-EARTH*r12i1p1*.nc', ##new
                     'MOHC_HadGEM2_regrid': cordex11_dir+'MOHC/*MOHC-HadGEM2-ES*r1i1p1*.nc', ##new

                     'RegCM-MPI_regrid': cordex11_dir + 'ICTP/*MPI*.nc', ## new
                     'RegCM_HadGEM2_regrid': cordex11_dir + 'ICTP/*HadGEM*.nc', ##new
                    #  'RegCM-MPI': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP/MPI-M-MPI-ESM-LR/ICTP-RegCM4-6/r1i1p1/pr*.nc',
                    #  'RegCM-HadGEM2': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP//MOHC-HadGEM2-ES/ICTP-RegCM4-6/r1i1p1/p*',
            }

    elif other == 'RegCM':
            runlist = {
                     # 'RegCM-MPI': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP/MPI-M-MPI-ESM-LR/ICTP-RegCM4-6/r1i1p1/pr*.nc',
                     # 'RegCM-HadGEM2': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX11/EUR-11/historical/day/pr/ICTP//MOHC-HadGEM2-ES/ICTP-RegCM4-6/r1i1p1/p*',
                     # 'RegCM-HadGEM_50km': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50//pr/ICTP/MOHC-HadGEM2-ES/ICTP-RegCM4-3/r1i1p1/*_v1_*.nc', ####
                     'RCA-MIROC_50km': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50/pr/SMHI/MIROC-MIROC5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                       }
 
    elif other == 'PRIMAVERA_orig':
        member_dict = {#'MPI-ESM1-2-HR': ['r1i1p1f1',], 
                   #'MPI-ESM1-2-XR': ['r1i1p1f1',],
                   #'HadGEM3-GC31-HM': ['r1i1p1f1', ],#'r1i2p1f1'],
				   #'HadGEM3-GC31-LL': ['r1i1p1f1',],# 'r1i2p1f1', 'r1i3p1f1', 'r1i4p1f1', 'r1i5p1f1', 'r1i6p1f1', 'r1i7p1f1', 'r1i8p1f1'], 
                   #'CNRM-CM6-1-HR': ['r1i1p1f2',],
                   #'CNRM-CM6-1': ['r1i1p1f2',],
                   #'EC-Earth3P-HR': ['r1i1p2f1',], # one more member available on request
                   #'EC-Earth3P': ['r1i1p2f1',], # one more member available on request
                   #'ECMWF-IFS-LR': ['r6i1p1f1',], #['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1', 'r7i1p1f1', 'r8i1p1f1',],
                   #'ECMWF-IFS-HR': ['r6i1p1f1',], #['r2i1p1f1', 'r3i1p1f1', 'r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1'],
                   #'CMCC-CM2-HR4': ['r1i1p1f1', ], 
                   'CMCC-CM2-VHR4': ['r1i1p1f1', ],
                      }       
        institution_dict = {'MPI-ESM1-2-HR': 'MPI-M', 'MPI-ESM1-2-XR': 'MPI-M', 'HadGEM3-GC31-HM': 'MOHC', 'HadGEM3-GC31-LL': 'MOHC',
                            'CNRM-CM6-1-HR': 'CNRM-CERFACS', 'CNRM-CM6-1': 'CNRM-CERFACS', 'EC-Earth3P-HR': 'EC-Earth-Consortium',
                            'EC-Earth3P': 'EC-Earth-Consortium', 'ECMWF-IFS-LR': 'ECMWF', 'ECMWF-IFS-HR':'ECMWF'} 
        runlist = {}
        for model in member_dict:
            for member in member_dict[model]:
                ### ME: change pr to ta here
                version = subprocess.check_output(['/home/users/croberts01/bin/primdata', '--latest_version', '--table','day', '--var','pr', '--exp', 'hist-1950', '--model', model, '--member', member])[:-1]
                if 'CMCC' in model:
                    runlist['{}_{}'.format(model, member)] = '/gws/nopw/j04/primavera3/cache/sberthou/CMCC_daily/pr_{}_hist-1950_{}_gn_*.nc'.format(model, member)
                else:
                    runlist['{}_{}'.format(model, member)] = primavera_path.format(institution_dict[model], model, member, version, model, member)
        print runlist
    elif other == 'PRIMAVERA': #regridded on CORDEX
        prim_cordex_path = '/gws/nopw/j04/primavera3/cache/sberthou/pr_{}_{}_regridded_on_{}.nc'
        cordex_grid = 'EUROCORDEX'
        member_dict = {'MPI-ESM1-2-XR': ['r1i1p1f1',],
                       'HadGEM3-GC31-HM': ['r1i1p1f1', ],#'r1i2p1f1'],
                       'CNRM-CM6-1-HR': ['r1i1p1f2',],
                       'EC-Earth3P-HR': ['r1i1p2f1',], # one more member available on request
                       'ECMWF-IFS-HR': ['r1i1p1f1', ], #'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1'],
                       'CMCC-CM2-VHR4': ['r1i1p1f1',], 
                      }
        runlist = {}
        for model in member_dict:
            for member in member_dict[model]:
                runlist['{}_{}_{}'.format(model, member, cordex_grid)] = prim_cordex_path.format(model, member, cordex_grid)
        print runlist

    elif other == 'PRIMAVERA_reduced':
        prim_cordex_path = '/gws/nopw/j04/primavera3/cache/sberthou/pr_{}_{}_regridded_on_{}.nc'
        cordex_grid = 'EUROCORDEX'
        member_dict = {'MPI-ESM1-2-XR': ['r1i1p1f1',],
                       'HadGEM3-GC31-HM': ['r1i1p1f1', ],#'r1i2p1f1'],
                       'CNRM-CM6-1-HR': ['r1i1p1f2',],
                       'EC-Earth3P-HR': ['r1i1p2f1',], # one more member available on request
                      # 'ECMWF-IFS-HR': ['r1i1p1f1', ], #'r2i1p1f1', 'r3i1p1f1', 'r4i1p1f1', 'r5i1p1f1', 'r6i1p1f1'],
                      # 'CMCC-CM2-VHR4': ['r1i1p1f1',],
                      }
        runlist = {}
        for model in member_dict:
            for member in member_dict[model]:
                runlist['{}_{}_{}'.format(model, member, cordex_grid)] = prim_cordex_path.format(model, member, cordex_grid)
        print runlist
    
    elif other == 'CMIP5_reduced':
        # Note, this is only used for the plotting in interface_plot_european...
        # the processing was done on another machine at ETH by M. E. Demory
             runlist =  {# 'CanESM2': '', #not for BI,AL
                     'CNRM-CM5': '',
                     'CSIRO-Mk3-6-0': '',
                     'EC-EARTH': '',
                     #'EC-EARTH12': '',
                     'GFDL-ESM2G': '', #not for BI
                     #'GFDL-ESM2M': '', #not for BI
                     'HadGEM2-ES': '', #not for BI
                     'IPSL-CM5A-LR': '', #not for BI
                     #'IPSL-CM5A-MR': '', #not for BI
                     'MIROC5': '',
                     'MPI-ESM-LR': '',
                     #'MPI-ESM-LR2': '',
                     'NorESM1-M': '', #not for BI
                    } 
    elif other == 'CMIP5':
              runlist = { 'ACCESS1-01': '', #not for BI
                          'ACCESS1-31': '', #not for BI
                          #'BNU-ESM': '', #not for AL,BI
                          'CCSM41': '',
                          'CESM1-BGC': '',
                          'CESM1-CAM5': '',
                          'CESM1-FASTCHEM1': '',
                          #'CMCC-CESM': '', #not for AL, BI, FR
                          'CMCC-CM': '',
                          'CMCC-CMS': '',
                          'CNRM-CM5': '',
                          'CSIRO-Mk3-6-0': '',
                          #'CSIRO-Mk3L-1-21': '',
                          #'CanCM41': '', #not for AL,BI
                          #'CanESM2': '', #not for AL,BI
                          'EC-EARTH': '',
                          #'FGOALS-g2': '',  #not for AL,BI
                          #'GFDL-CM3': '',   #pb to concatenate
                          'GFDL-ESM2G': '', #not for BI
                          'GFDL-ESM2M': '', #not for BI
                          'GISS-E2-H2': '', #not for BI
                          'GISS-E2-R': '',  #not for BI
                          'HadCM3': '',     #not for BI,MD,NEE
                          'HadGEM2-AO': '', #not for BI
                          'HadGEM2-CC': '', #not for BI
                          'HadGEM2-ES': '', #not for BI
                          'IPSL-CM5A-LR': '', #not for BI
                          'IPSL-CM5A-MR': '', #not for BI
                          'IPSL-CM5B-LR': '', #not for BI
                         # 'MIROC-ESM': '',    #not for AL,BI
                         # 'MIROC-ESM-CHEM': '', #not for AL,BI
                          'MIROC4h': '',
                          'MIROC5': '',
                          'MPI-ESM-LR': '',
                          'MPI-ESM-MR': '',
                          'MPI-ESM-P': '',
                          'MRI-CGCM3': '',
                          'MRI-ESM1': '',
                          'NorESM1-M': '',   #not for BI
                          #'bcc-csm1-1': '',  #not for AL,BI
                          'bcc-csm1-1-m': '',
                          'inmcm4': '',      #not for BI
                 }
    else:
        runlist = {
              'eobs_cordex50':'/gws/nopw/j04/primavera3/cache/sberthou/E-OBS/rr_0.50deg_reg_v17.0_EUROCORDEX.nc'
              #'eobs': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/E-OBS/rr_0.50deg_reg_v17.0.nc',
              #country+'_cordex50':runlist_frequency[frequency][country+'_cordex50'],
        }
   
    return runlist, maskedreglist
