# import ucc.contrib.fcfa.plotting as fcplt
import subprocess


class ParaArgs(object):
    '''
    Defines a class to collate the calling arguments to functions run in parallel with the
    pool method.
    Usage: Loop as in serial processing to generate a list of all calling arguments.
           outsiede the loop, open a parallel pool region and map the list of calling
           arguments
    '''

    def __init__(self, cmd=None, msg=None):
        self.cmd = cmd
        self.msg = msg


def shellcmd(args):
    try:
        retcode = subprocess.call(args.cmd, shell=True)
        if retcode < 0:
            print('syst.cmd terminated by signal', retcode)
        elif retcode:
            print('syst.cmd returned in ', args.msg, '', retcode)
    except OSError as ex:
        print("Execution failed in " + args.msg + ": ", ex)

##### these are the name of the regions defined in dict_subregions in european_masked_subregion.py ####
##### all the subregions defined in this dictionary entry will be run ###########

countries = [  # 'france', 'alps',
    'prudence', 
#     'uk',
#    'alps', 'ibp', 'ca',
#    'france',
    # 'germany', 'switzerland', 'spain',
] 
frequencies = ['d', ] # daily
resolutions = ['False', ]  # don't care for CORDEX analysis - this was used in another version - leave it so.
bins = ['exponential100',] # name of the bins defined in process_data_europe.py # exponential100 is used in the GMD paper
                           # other options are traditional, martin, exponential, martin_100. Martin_100 is the one used by Klingaman et al. (2018)
wet_dry_spells = None  # wet or dry # this is not used for the GMD paper, leave as "None"
thrs_spells = 0.1 # only used when wet_dry_spells is "wet" or "dry"

for country in countries:
    for freq in frequencies:
        for res in resolutions:
            for bin1 in bins:
                #  make sure the right python is called here: python 2.?
                # this loop calls another python script which launches jobs on the cluster
                cmd = 'python2.7 launch_european_masked_subregion.py -c {} -r {} -f {} -s djf mam jja son -b {} -e True -w {} -t {} -o new_CORDEX_11_regrid'.format( #-o PRIMAVERA_regridded_CORDEX
                    country,
                    res,
                    freq,
                    bin1,
                    wet_dry_spells,
                    thrs_spells)
                argument = ParaArgs(cmd=cmd,
                                    msg='{}_{}_{}_{} failed'.format(country, res, freq, bin1))
                shellcmd(argument)
