# import ucc.contrib.fcfa.plotting as fcplt
import subprocess

countries = [
    'prudence',]
#     'france',]# 'alps', 'nee', 'ca', 'ibp', 'uk', 'sc', 'ce', 'md']#['mali','niger', 'benin']#['nimrod', 'germany', 'switzerland']#, 'niger', 'benin']
frequencies = ['d']
resolutions = ['False']  # ['True=n512', 'False=mi-ao438'] ## M-E, don't worry about this
bins = ['exponential100']  # ['martin']
seasons = ('djf', 'mam', 'jja', 'son')  # ('amjjas',)#''djf', 'mam', 'jja', 'son')


def shellcmd(cmd, msg_func):
    try:
        retcode=subprocess.call(cmd, shell=True)
        if retcode<0:
            print('syst. cmd terminated by signal', -retcode)
        elif retcode:
            print('syst. cmd returned in ',msg_func,' ', retcode)
    except OSError as ex:
        print("Execution failed in "+msg_func,+": ",ex)


for country in countries:
    for freq in frequencies:
        for res in resolutions:
            for bin1 in bins:
#                for season in seasons:
                    ### check python interpreter is right (python 2.?)

                    # important: -n can be: contrib, frac, normed or normed_wet
                    # contrib=actual contribution (bin_mean*frequency)
                    # frac= contrib / mean precip over the region
                    # normed=frequency
                    # normed_wet=frequency only over wet days (1mm.day threshold)

                    # note the -w wet and -t 0.1 options are ignored if -g is set to 1d. if set to 2d, then they are taken into account

                    cmd = 'python2.7 interface_plot_european_masked_subregions.py -c {} -r {} -f {} -s {} -b {} -n contrib -g 1d -w wet -t 0.1'.format(
                        country,
                        res,
                        freq,
                        ' '.join([s for s in seasons]),
                        bin1)
                    shellcmd(cmd, '{}_{}_{}_{} failed'.format(country, res, freq, bin1))
