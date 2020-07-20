
import european_masked_subregion as ems
import __init__
import iris 
import iris.analysis as iana
import re
from matplotlib.colors import from_levels_and_colors
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.cm as cm
import cartopy.crs as ccrs
import numpy as np
from scipy import stats
from statsmodels.sandbox.stats import multicomp

def main():
    datadir = '/home/users/sberthou/PrecipExtremes/'
    frequency = 'd'
    full_ensemble = False
    if full_ensemble:
     runlist, _ = ems.get_runlist_region(frequency, False, 'prudence', 'PRIMAVERA')
     PRIMAVERA_ens = runlist.keys()
     runlist1, _ = ems.get_runlist_region(frequency, False, 'prudence', 'CORDEX-44')
     CORDEX_ens = runlist1.keys()
     runlist2, _ = ems.get_runlist_region(frequency, False, 'prudence', 'CORDEX11')
     CORDEX_ens_11 = runlist2.keys()
     runlist3, _ = ems.get_runlist_region(frequency, False, 'prudence', 'CORDEX-11_regrid')
     CORDEX_ens11_regrid = runlist3.keys()
    else:
     runlist, _ = ems.get_runlist_region(frequency, False, 'prudence', 'PRIMAVERA_reduced')
     PRIMAVERA_ens = runlist.keys()
     runlist1, _ = ems.get_runlist_region(frequency, False, 'prudence', 'CORDEX-44_reduced')
     CORDEX_ens = runlist1.keys()
     runlist3, _ = ems.get_runlist_region(frequency, False, 'prudence', 'CORDEX-11_regrid_reduced')
     CORDEX_ens11_regrid = runlist3.keys()
 
    obs_ens = ['eobs_cordex50', 'ca_cordex50', 'france_cordex50', 'ibp_cordex50', 'alps_cordex50', 'uk_cordex50',]
    seasons = ['jja', 'djf', 'mam', 'son']
    ref_grid = iris.load_cube('{}/mean_nc/ibp_cordex50_jja_mean_2000.nc'.format(datadir))
    ensembles = {'OBS':obs_ens, 'PRIMAVERA': PRIMAVERA_ens, 
                 'EUR-44': CORDEX_ens, #'CORDEX-11': CORDEX_ens_11, 
                 'EUR-11': CORDEX_ens11_regrid}
    for season in seasons:
        print season
        cubelist_dict = {}
        for ensemble in ensembles:
            print ensemble
            cubelist = iris.cube.CubeList([])
            for num, simu in enumerate(ensembles[ensemble]):
                print simu
                cubelist0 = iris.load('{}/mean_nc/{}_{}_mean*.nc'.format(datadir, simu, season))
                if len(cubelist0) == 1:
                    cube = cubelist0[0]
                else:
                    cube = cubelist0.extract('precipitation_flux')[0]
                modify_grid(cube, simu)
                coord_names = [coord.name() for coord in cube.dim_coords]
                print cube
                try:
                    lonname, = filter(lambda x: re.search(r'[Ll]on|x', x), coord_names)
                except:
                    __init__.add_CORDEX_grid(cube, 'EUROCORDEX')
                
                newcube = cube.regrid(ref_grid, iana.Linear())
                if ensemble is not 'OBS':
                   member_coord = iris.coords.AuxCoord([num],
                                              long_name='member_number',
                                              units='1')
                   newcube.add_aux_coord(member_coord)
                   newcube.attributes = None
                   newcube.cell_methods = None
                   try:
                     newcube.remove_coord('time')
                   except:
                     pass
                newcube1 = precip_to_mm(newcube)
                cubelist.append(newcube1)
            cubelist_dict[ensemble] = cubelist.merge()
            print cubelist_dict[ensemble]
        obs = combine_obs(cubelist_dict['OBS'])
        add_reduced = ''
        if not full_ensemble:
           add_reduced = '_reduced'
        figfile = '{}/images/mean_bias{}_{}'.format(datadir, add_reduced, season)
        plot_maps(obs, cubelist_dict, plot_order=['PRIMAVERA', 'EUR-11', 'EUR-44'], figfile=figfile, season=season)
       # save_ensemble_mean(obs, cubelist_dict, season, datadir, add_reduced)


def combine_obs(obs_cubelist):
    obs = obs_cubelist[0].copy()
    obs0 = obs.copy()
    for num, new_obs in enumerate(obs_cubelist[1:]):
        #new_obs1 = new_obs.collapsed('time', iana.MEAN) #mean over all the years
        obs.data[~new_obs.data.mask] = new_obs.data[~new_obs.data.mask]
#    plt.figure()
#    qplt.pcolormesh(obs-obs0, cmap=cm.get_cmap('RdBu_r'), vmin = -1, vmax=1)
#    plt.savefig('obsdiff.png')
    return obs


def calculate_sig_val(obs, ensemble):
    sig_val = obs[0,:, :].copy()
    _, pval = stats.ttest_ind(obs.data, ensemble.data, axis=0, equal_var=False)
    print pval.shape
    print sig_val.shape
    print pval
    sig_val.data = pval
    sig_val = mask_cube_where(sig_val, obs[0, :, :].data.mask)
    return sig_val


def mask_cube_where(cube, cdt):
        '''
        Returns a cube with masked data where the condition cdt is met

        :param iris.cube cube:
        :param bool cdt: can be a condition on 3d data of the cube
        :return iris.cube masked_cube: masked cube
        '''
        masked_cube = cube.copy()
        if not np.ma.is_masked(masked_cube.data,):
        # if not isinstance(masked_cube.data, np.ma.core.MaskedArray):
            masked_cube.data = np.ma.asarray(cube.data)  # np.ma.MaskedArray(cube.data,  mask = False)
        masked_cube.data = np.ma.masked_where(cdt, masked_cube.data, copy=True)
        return masked_cube


def plot_maps(obs, cubelist_dict, plot_order, figfile, season):
    fig = plt.figure(figsize=(3.5, 3))
    axx = fig.add_subplot(221, projection=ccrs.PlateCarree())
    plot1 = iplt.pcolormesh(obs.collapsed('time', iana.MEAN), cmap=cm.get_cmap('terrain_r'), vmin=0, vmax=7)
    plt.title('a. observations ({})'.format(season.upper()), fontsize=8)
    cbar = plt.colorbar(plot1, cax=plt.gcf().add_axes((0.02, 0.1, 0.43, 0.02)),
                                        orientation='horizontal',)
    letter_dict = {'PRIMAVERA': 'b',
                   'EUR-44': 'd',
                   'EUR-11': 'c'}
    cbar.ax.tick_params(labelsize=8)
    plt.title('mm/day', fontsize=8, y=0.7)
    axx.set_extent([-22, 41, 33, 70])
    for num, dataset in enumerate(plot_order, 2):
       numnum = '22'+str(num)
       axx = fig.add_subplot(int(numnum), projection=ccrs.PlateCarree())
       model = cubelist_dict[dataset][0].collapsed('member_number', iana.MEAN)
       #plot = iplt.pcolormesh(model-obs, cmap=cm.get_cmap('RdBu'), vmin = -2, vmax=2)
       levels = np.linspace(-2, 2, 11)
       extend = 'both'
       sig_val = calculate_sig_val(obs, cubelist_dict[dataset][0])
       cube2plot = model-obs.collapsed('time', iana.MEAN)
       cmap=cm.get_cmap('RdBu', 12)
       colors = np.concatenate([cmap(range(0, 12)), ])
       cmap1, norm = from_levels_and_colors(levels, colors, extend='both')
       plot = iplt.pcolormesh(cube2plot, norm=norm, cmap=cmap1)
       global_mask = generate_global_pvalue_mask(sig_val.data, sig_lev=0.1)
       cube2plot = mask_cube_where(cube2plot, global_mask)
       iplt.contourf(cube2plot, levels=levels, cmap=cm.get_cmap('RdBu'),
                      extend=extend, hatches=['.....'], zorder=10,
                      alpha=0.5)

       plt.title('{}. {} - observations'.format(letter_dict[dataset], dataset), fontsize=8)
       axx.set_extent([-22, 41, 33, 70])
    cbar = plt.colorbar(plot, cax=plt.gcf().add_axes((0.52, 0.1, 0.43, 0.02)),
                                        orientation='horizontal',)
    cbar.ax.tick_params(labelsize=8)
    plt.title('mm/day', fontsize=8, y=0.7)
    plt.box(False)
    plt.tight_layout(rect=[0, 0.12, 1, 1], pad=0.05)
    plt.savefig(figfile+'_'.join(r for r in plot_order) +'.png', dpi=300)


def save_ensemble_mean(obs, cubelist_dict, season, datadir, add_reduced):
    for dataset in cubelist_dict:
        if dataset is not 'OBS':
            model = cubelist_dict[dataset][0].collapsed('member_number', iana.MEAN)
            iris.save(model, '{}/ens_mean_nc/{}{}_{}.nc'.format(datadir, dataset, add_reduced, season))
    #iris.save(obs, '{}/ens_mean_nc/{}_{}.nc'.format(datadir, 'obs', season))


def modify_grid(r_ppn, simu):
    #if ('ALADIN' in simu) | ('ALARO' in simu):
    #      r_ppn.coord('projection_x_coordinate').convert_units('m')
    #      r_ppn.coord('projection_y_coordinate').convert_units('m')
    if ('WRF-EC-EARTH_50km' in simu) | ('WRF-IPSL_50km' in simu):
          r_ppn.remove_coord('grid_latitude')
          r_ppn.remove_coord('grid_longitude')
          __init__.add_CORDEX_grid(r_ppn)
    if 'RegCM' in simu:
       if '50km' in simu:
          cube_ref = iris.load_cube('/group_workspaces/jasmin4/upscale/cache/demory/Segolene/CORDEX50//pr/ICTP/MOHC-HadGEM2-ES/ICTP-RegCM4-3/r1i1p1/pr_EUR-44_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-3_v1_day_1950010112-1950123012.nc')
          el_axis = cube_ref.coord('projection_x_coordinate').coord_system.ellipsoid.semi_major_axis
          cube_ref.coord('projection_x_coordinate').coord_system.ellipsoid.semi_minor_axis = el_axis
          cube_ref.coord('projection_y_coordinate').coord_system.ellipsoid.semi_minor_axis = el_axis
          r_ppn.add_dim_coord(cube_ref.coord('projection_x_coordinate'), 1)
          r_ppn.add_dim_coord(cube_ref.coord('projection_y_coordinate'), 0)


def precip_to_mm(cube, freq='day'):
    '''
    What it does  : Convert Precip Units to mm/day
    How to use it : precip_to_mm(cube)
    Example       : cube = precip_to_mm(cube)
    :param str freq: day or 3h for example  
    '''

    if ( cube.units in ['kg m-2 s-1', 'kg m-2'] ):
        rho = iris.coords.AuxCoord(1000.0,long_name='ref density',units='kg m-3')
        outcube = cube/rho
        if cube.units == 'kg m-2 s-1':
            outcube.convert_units('mm/{}'.format(freq))
        else:
            outcube.convert_units('mm')
            outcube.units = 'mm/{}'.format(freq)
        outcube.rename(cube.name())
    
    elif ( cube.units == 'm s-1' ):
        outcube = cube.copy()
        outcube.convert_units('mm/{}'.format(freq))
        outcube.standard_name = cube.standard_name
    elif ( cube.units in ['mm 3h-1', 'mm (3h)-1'] ):
        outcube = cube.copy()
        outcube.convert_units('mm/{}'.format(freq))
        outcube.standard_name = cube.standard_name
        warnings.warn('Converted units from mm (3h)-1 to {}'.format(outcube.units))
    elif ( cube.units in ['mm {}-1'.format(freq), 'mm', 'mm.d-1', 'mm ({})-1'.format(freq)] ):
        outcube = cube.copy()
        outcube.units = 'mm/{}'.format(freq)
    elif ( cube.units == 'mm/90days' ):
        outcube = cube.copy()
    elif ( cube.units == 'mm (3h)-1', 'mm 3h-1' ):
        outcube = cube.copy()
        outcube.convert_units('mm/{}'.format(freq))
        outcube.standard_name = cube.standard_name
    else:
        raise ValueError('cube not in mm/{}, kg/m2/s or m/s '.format(freq), cube.units
                         )
    outcube.long_name = None
    outcube.standard_name = 'precipitation_flux'
    return outcube


def generate_global_pvalue_mask(p_val, sig_lev=0.05):
     '''
     This function calculates the local pvalues from the bootstrapping
     (calling estimate_local_pvalue_from_bootstrap), and then calculates the global p_value with
     multicomp.multipletests and generates a mask of significant values
     (then can be applied to the diff cube and plotted using hatches = ["..."] )

     :param iris.cube control_cube: 2D cube
     :param iris.cube control_cube_simulated: 3D cube with 'simulated_sample' coordinate
     :param iris.cube future_cube: 2D cube
     :param iris.cube future_cube_simulated: 3D cube with 'simulated_sample' coordinate#
     :param float sig_lev: global level of significance wanted to be achieved

     :return annual_sig_test_mask
     '''
     alpha_fdr_adj = 2.0

     p_val_compressed = p_val.data.compressed() if isinstance(p_val.data, np.ma.MaskedArray) else p_val.data.flatten()
     p_val2 = multicomp.multipletests(p_val_compressed,
                                      alpha=alpha_fdr_adj * sig_lev / 2.0,
                                      method='fdr_bh', is_sorted=False,
                                      returnsorted=False)
     annual_sig_test_mask_raw = np.logical_not(p_val2[0])

     if isinstance(p_val.data, np.ma.MaskedArray):  # reshape and remask
         annual_sig_test_mask = np.empty_like(p_val.data)
         np.place(annual_sig_test_mask, annual_sig_test_mask.mask, 1)
         np.place(annual_sig_test_mask, ~annual_sig_test_mask.mask,
                  annual_sig_test_mask_raw)
     else:
         annual_sig_test_mask = np.reshape(annual_sig_test_mask_raw, p_val.data.shape)

#     print(p_val2[0][0:10], p_val2[1][0:10], p_val_compressed[0:10])
#     print(p_val.data.max(), p_val.data.min(), p_val.data.shape)
#     print(p_val2[1].max(), p_val2[1].min(), np.median(p_val2[1]), np.mean(p_val2[1]), p_val2[1].shape)

     return annual_sig_test_mask


if __name__ == '__main__':
    main()        
