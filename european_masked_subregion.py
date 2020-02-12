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
                                    ('NEE', datadir+'mask_PRUDENCEreg_NEE.nc'),
                                    #('IP', datadir+'mask_PRUDENCEreg_IP.nc'), # Once BI region works, then uncomment all the others here
                                    #('FR', datadir+'mask_PRUDENCEreg_FR.nc'),
                                    ('MD', datadir+'mask_PRUDENCEreg_MD.nc'),
                                    #('CA', datadir+'mask_PRUDENCEreg_CA.nc'),
                                    ('CE', datadir+'mask_PRUDENCEreg_CE.nc'),
                                    ('SC', datadir+'mask_PRUDENCEreg_SC.nc'),
                                    #('AL', datadir+'mask_PRUDENCEreg_AL.nc'),
                                    ##'EE': datadir+'mask_PRUDENCEreg_EE.nc',
                                    #('BI', datadir+'mask_PRUDENCEreg_BI.nc'),
                                    #'CH': datadir+'mask_PRUDENCEreg_CH.nc',
                                    # 'TCZ': datadir+'mask_PRUDENCEreg_TCZ.nc'
                                    ]),
                        'uk': {'BI': datadir+'mask_PRUDENCEreg_BI.nc',},
                        'ca': {'CA': datadir+'mask_PRUDENCEreg_CA.nc',},
                        'alps': {'AL': datadir+'mask_PRUDENCEreg_AL.nc',},
                        'ibp': {'IP': datadir+'mask_PRUDENCEreg_IP.nc',},
                        'france': {'FR': datadir+'mask_PRUDENCEreg_FR.nc',},
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

    elif other == 'CORDEX':
        cordex_dir = '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50/'
        runlist = {
             #EUR44 models
                     #'CCLM-MPI_50km': cordex_dir+'/pr/CLMcom/MPI-M-MPI-ESM-LR/CLMcom-CCLM4-8-17/r1i1p1/*_v1_*.nc',
                     #'CCLM5-CNRM_50km': cordex_dir+'/pr/ETH/CNRM-CERFACS-CNRM-CM5/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     #'CCLM5-EC-EARTH_50km': cordex_dir+'/pr/ETH/ICHEC-EC-EARTH/CLMcom-CCLM5-0-6/r12i1p1/*_v1_*.nc',
                     #'CCLM5-HadGEM_50km': cordex_dir+'/pr/ETH/MOHC-HadGEM2-ES/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc', #problem? -> check
                     #'CCLM5-MIROC_50km': cordex_dir+'/pr/ETH/MIROC-MIROC5/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     #'CCLM5-MPI_50km': cordex_dir+'/pr/ETH/MPI-M-MPI-ESM-LR/CLMcom-CCLM5-0-6/r1i1p1/*_v1_*.nc',
                     #'HIRHAM-EC-EARTH_50km': cordex_dir+'/pr/DMI/ICHEC-EC-EARTH/DMI-HIRHAM5/r3i1p1/*_v1_*.nc',
                     #'WRF-IPSL_50km': cordex_dir+'/pr/IPSL-INERIS/IPSL-IPSL-CM5A-MR/IPSL-INERIS-WRF331F/r1i1p1/*_v1_*.nc', #problems
                     #'WRF-EC-EARTH_50km': cordex_dir+'/pr/NUIM/ICHEC-EC-EARTH/NUIM-WRF341E/r1i1p1/*_v1_*.nc', #problems!!!!check it worked
                     #'RACMO-EC-EARTH_50km': cordex_dir+'/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r1i1p1/*_v1_*.nc',
                     #'RACMO-EC-EARTH12_50km': cordex_dir+'/pr/KNMI/ICHEC-EC-EARTH/KNMI-RACMO22E/r12i1p1/*_v1_*.nc', ######
                     #'RACMO-HadGEM_50km': cordex_dir+'/pr/KNMI/MOHC-HadGEM2-ES/KNMI-RACMO22E/r1i1p1/*_v2_*.nc',
                     #'REMO-MPI_50km': cordex_dir+'/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r1i1p1/*_v1_*.nc',
                     #'REMO-MPI2_50km': cordex_dir+'/pr/MPI-CSC/MPI-M-MPI-ESM-LR/MPI-CSC-REMO2009/r2i1p1/*_v1_*.nc',#####
                     #'RCA-CanESM_50km': cordex_dir+'/pr/SMHI/CCCma-CanESM2/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-CNRM_50km': cordex_dir+'/pr/SMHI/CNRM-CERFACS-CNRM-CM5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-CSIRO_50km': cordex_dir+'/pr/SMHI/CSIRO-QCCCE-CSIRO-Mk3-6-0/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-EC-EARTH_50km': cordex_dir+'/pr/SMHI/ICHEC-EC-EARTH/SMHI-RCA4/r12i1p1/*_v1_*.nc',
                     #'RCA-IPSL_50km': cordex_dir+'/pr/SMHI/IPSL-IPSL-CM5A-MR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-MIROC_50km': cordex_dir+'/pr/SMHI/MIROC-MIROC5/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-HadGEM_50km': cordex_dir+'/pr/SMHI/MOHC-HadGEM2-ES/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-MPI_50km': cordex_dir+'/pr/SMHI/MPI-M-MPI-ESM-LR/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-NorESM_50km': cordex_dir+'/pr/SMHI/NCC-NorESM1-M/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'RCA-GFDL_50km': cordex_dir+'/pr/SMHI/NOAA-GFDL-GFDL-ESM2M/SMHI-RCA4/r1i1p1/*_v1_*.nc',
                     #'ALADIN-CNRM1_50km': cordex_dir+'/pr/CNRM/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN53/r1i1p1/*_v1_*.nc',
                     #'ALADIN-CNRM2_50km': cordex_dir+'/pr/HMS/CNRM-CERFACS-CNRM-CM5/HMS-ALADIN52/r1i1p1/*_v1_*.nc',
                     'WRF-CanESM_50km': cordex_dir+'/pr/UCAN/CCCma-CanESM2/UCAN-WRF341I/r1i1p1/*_v2_*.nc', ####
                     #'RegCM-HadGEM_50km': cordex_dir+'/pr/ICTP/MOHC-HadGEM2-ES/ICTP-RegCM4-3/r1i1p1/*_v1_*.nc', ####
                     #'ALARO-CNRM_50km': cordex_dir+'/pr/RMIB-UGent/CNRM-CERFACS-CNRM-CM5/RMIB-UGent-ALARO-0/r1i1p1/*_v1_*.nc',
                   }
    elif other == 'PRIMAVERA':
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
    elif other == 'PRIMAVERA_regridded_CORDEX':
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

    else:
        runlist = {
              'eobs_cordex50':'/gws/nopw/j04/primavera3/cache/sberthou/E-OBS/rr_0.50deg_reg_v17.0_EUROCORDEX.nc'
              #'eobs': '/group_workspaces/jasmin4/upscale/cache/demory/Segolene/E-OBS/rr_0.50deg_reg_v17.0.nc',
            #  country+'_cordex50':runlist_frequency[frequency][country+'_cordex50'],
        }

    return runlist, maskedreglist
