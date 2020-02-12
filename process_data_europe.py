'''
Created on Jun 29, 2016

@author: cbates

Generating histogram data from given model data
'''
import __init__
import os
import argparse
import sys
import iris
import datetime
import numpy as np
import subprocess

iris.FUTURE.cell_datetime_objects = True

### change datadir to where you want files to be output: they will be in datadir/jsondir
config = {'datadir': "/home/users/sberthou/PrecipExtremes/", #/data/local/hadmm/upscale_precip/daily/",
    'jsondir': "srex_pdf_json/",
    'lsmasks': {} # not used here
}

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    # User input of histfiles and infiles
    parser.add_argument("-o", type=str, help="Last part of histogram file")
    parser.add_argument("-i", nargs="*", type=str, help="Input file ")
    # User input of desired region, if not specified then generate data for all regions
    parser.add_argument("-r", nargs="*", type=str, help="Region ", default=__init__.srex_reg.keys())
    parser.add_argument("-n", type=str, help="Resolution")
    parser.add_argument("-g", type=str, help="Regrid resolution", default=None)
    parser.add_argument("-s", nargs='*', type=str, help="Season list", default=[None])
    parser.add_argument("-q", type=str, help="Season", default=None)
    parser.add_argument("-d", type=str, help="traditional or martin distribution?", default='traditional')
    parser.add_argument("-t", type=float, help="threshold for spell calculation (in mm/day)", default=0)
    parser.add_argument("-w", type=str, help="wet or dry spell calculation (wet/dry)", default=None)
    parser.add_argument("-m", type=str, help="mask file location", default=None)
    parser.add_argument("-a", type=str, help="add rain and snow for euro2.2km", default=None)
    parser.add_argument("-f", type=str, help="frequency (d, h, 3h, 15mn)", default='d')

    args = parser.parse_args()

    #    if len(args.i) == 1 and not os.path.exists(args.i[0]):
    #        print "\nNo files found to process\n"
    #        sys.exit(255)

    if args.i is None or args.o is None or args.r is None or args.n is None or args.s is None or args.d is None:
        parser.print_help()
        print "\nPlease provide valid arguments\n"
        exit()

    infiles = args.i
    histfile = args.o
    region = [args.r, ][0]
    resolution = [args.n, ]
    regrid_res = args.g
    season_list = args.s
    season = args.q
    distribution = args.d
    mask2d = args.m
    rain_and_snow = args.a
    config['infiles'] = infiles
    config['histfile'] = histfile
    thrs = args.t
    wet_or_dry_spells = args.w
    if wet_or_dry_spells == 'None':
        wet_or_dry_spells = None
    frequency = args.f

    dict_freq = {'h': 24.,
                 '3h': 8.,
                 'd': 1.,
                 '15mn': 24. * 4}

    jsondir_freq = {'d': 'srex_pdf_json/',
                    'h': 'srex_pdf_json_hourly/',
                    '3h': 'srex_pdf_json_3hourly/',
                    '15mn': 'srex_pdf_json_15mn/'}

    config['jsondir'] = jsondir_freq[frequency]

    ### I set the time constraint to 1950-2005, but this is where you can change it if needs be.
    time_con = iris.Constraint(time=lambda cell: 1971 <= cell.point.year < 2005)
    hourly2daily = False
    duration_bins = None

    #### this is where the bin definition is set
    if distribution == 'traditional':
        bins1 = np.arange(1500)
        bins = np.pad(bins1, (0, 1), 'constant', constant_values=3000) / 86400.
    elif distribution == 'martin':
        bin2 = np.exp(
            np.log(0.005) + np.sqrt(np.linspace(0, 99, 200) * ((np.square(np.log(120.) - np.log(0.005))) / 59.)))
        bins = np.pad(bin2, (1, 0), 'constant', constant_values=0) / 86400.
    elif distribution == 'martin_100':
        bin2 = np.exp(
            np.log(0.005) + np.sqrt(np.linspace(0, 99, 100) * ((np.square(np.log(120.) - np.log(0.005))) / 59.)))
        bins = np.pad(bin2, (1, 0), 'constant', constant_values=0) / 86400.
    elif distribution == 'exponential':
        bin2 = np.exp(np.log(0.02) + 0.12*np.linspace(0, 99, 200))
        bins = np.pad(bin2, (1, 0), 'constant', constant_values=0) / 86400.
    elif distribution == 'exponential100':
        bin2 = np.exp(np.log(0.02) + 0.12*np.linspace(0, 99, 100))
        bins = np.pad(bin2, (1, 0), 'constant', constant_values=0) / 86400.
    elif distribution == 'exponential40':
        bin2 = np.exp(np.log(0.2) + 0.4*np.linspace(0, 22, 40))
        bins = np.pad(bin2, (1, 0), 'constant', constant_values=0) / 86400.


    if wet_or_dry_spells == 'wet':
        #        bins1 = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0])
        #        bins = np.pad(bins1, (0,1), 'constant', constant_values = 5000)/86400.
        dlim = 60 * dict_freq[frequency]  # days, duration max
        dbins = np.arange(1, dlim)
        duration_bins = np.pad(dbins, (0, 1), 'constant', constant_values=int(1000 * dict_freq[frequency]))
    else:
        # dry spells
        dlim = 360 * dict_freq[frequency]  # days, duration max
        ilim = 0.99  # mm/day, intensity max
        dbins = np.arange(1, dlim)
        duration_bins = np.pad(dbins, (0, 1), 'constant', constant_values=1000)
        ibins = np.arange(0, ilim, 0.02)
        bins = np.pad(ibins, (0, 1), 'constant', constant_values=1)
        bins /= 86400.


    PreExt = __init__.PrecipExtremes(__init__.srex_reg, config, resolution, bins, duration_bins)
    for res in PreExt.resolutions:
        print 'loading, ', res
        PreExt.load_mask_ppn(res, maskfile=mask2d, rain_and_snow=rain_and_snow, time_con=time_con,
                             season=season, season_list=season_list)
        for r in region:
            print r
            # for season in season_list:
            print season
            targdir = '{}{}'.format(config['datadir'], config['jsondir'])
            if not os.path.isdir(targdir):
                cmd = 'mkdir -p {}'.format(targdir)
                shellcmd(cmd, 'failed to create {}'.format(targdir))
            targfile = "{}/{}_{}_{}_{}".format(targdir, res, r, season, config['histfile'])
            if not os.path.exists(targfile):
                #### histogram calculation
                PreExt.gendata(regrid_res, r, wet_or_dry_spells, thrs, hr2day=hourly2daily)
                #### histogram saving into json files in json_dir
                PreExt.dump_data(targfile)
            else:
                print "file {} exists. skipping".format(targfile)


def shellcmd(cmd, msg_func):
    try:
        retcode=subprocess.call(cmd,shell=True)
        if retcode<0:
            print('syst. cmd terminated by signal', -retcode)
        elif retcode:
            print('syst. cmd returned in ',msg_func,' ', retcode)
    except OSError as ex:
        print("Execution failed in "+msg_func,+": ",ex)

