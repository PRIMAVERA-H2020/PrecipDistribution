'''
Created on Sep 13, 2016

This routine launches process_data_europe.py jobs on the cluster.
It is called by launch_all_europe.py
One job per bintype, season and subregion

@author: sberthou
'''
import os
#import my_library
# import multiprocessing as mp
import european_masked_subregion as ems
import argparse
import subprocess

parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", type=str, help="for which country? (Alps/France/Spain/ukcpobs/nimrod/mali/niger/benin")
parser.add_argument("-r", type=str, help = "for the files regridded on n512 (or model for ammacatch ", default = 'False')
parser.add_argument("-f", type=str, help = "", default = 'd')
parser.add_argument("-o", type=str, help = "Other (ukcpm_runs, ammacatch", default = '')
parser.add_argument("-s", nargs='*', type=str, help = "Season list", default = [None])
parser.add_argument("-b", nargs='*', type = str, help = "traditional or/and martin distribution?", default = ['traditional'])
parser.add_argument("-e", type=str, help = "run on spice: True/False", default = 'False')
parser.add_argument("-w", type=str, help = "compute 2d hists: wet/dry", default = None)
parser.add_argument("-t", type=str, help = "threshold for calculating wet/dry spells (mm)", default = 0)

args = parser.parse_args()

if args.c is None or args.r is None or args.f is None or args.o is None or args.s is None or args.b is None or args.e is None:
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
season_list = args.s
print type(season_list), season_list
bintypes = args.b
run_on_spice = bool_from_string(args.e)
wet_dry_2dhist = args.w
if wet_dry_2dhist == 'None':
   wet_dry_2dhist=None
thrs=args.t

# get the runlist and the region masks from european_masked_subregion.py
runlist, maskedreglist = ems.get_runlist_region(frequency, n512, country, other=other)

# batch job characteristics:
spice_freq = {'mem': {'h': '40', 'd': '20'}, ### ME: this is where you can change how much memory you use, 20GB by default
              'queue': {'h': 'normal', 'd': 'normal'},
              'time': {'h': '180', 'd': '100'},  ### ME: this is where you can change time limit, default is 100seconds
}

class para_args(object):
    '''
    Defines a class to collate the calling arguments to functions run in parallel with the 
    pool method.
    Usage: Loop as in serial processing to generate a list of all calling arguments.
           outsiede the loop, open a parallel pool region and map the list of calling 
           arguments
    Example: import FCFA_Utils_SpeedUp_and_SaveMemory as f_para
             import multiprocessing as mp
             for scenario is scen:
                 for season in seasons:
                     for model in models:
                         call_args.append(f_para.hycristal_args(scenario= scenario,
                                          season=season, model=models, out_dir = out_dir)
             p = mp.Pool(n_procs)
             p.map(parallel_function,call_args)
             
             def parallel_function(ca):
                  result_field = do_whatever(f_para.ca.model,f_para.ca.model,f_para.ca.season)
                  out_filename = create_filename(f_para.ca.model,f_para.ca.model,f_para.ca.season)
                  iris.fileformats.netcdf.save(result_field, out_dir + outfilename)

    Notes: In all the processes using this class, the parallel region ends by writing into a 
           file, so returning values to the serial region is not tested.
           
           I reckon that when every parallel region ends, it frees the memory it has used, 
           but if you know better, please let me know.
                                          
    Who : Jorge Bornemann
    When : 12/09/2016
    '''

    def __init__(self, cmd=None, msg=None):
        self.cmd = cmd
        self.msg = msg


def shellcmd(args):
    try:
            retcode = subprocess.call(args.cmd, shell=True)
            if retcode < 0:
                    print 'syst.cmd terminated by signal', retcode
            elif retcode:
                    print 'syst.cmd returned in ',args.msg,'',retcode
    except OSError as ex:
            print "Execution failed in "+args.msg+": ",ex


call_args=[]
for bintype in bintypes:
   print bintype
   for reg in maskedreglist.keys():
     print reg
     for s in season_list:
       for res in runlist.keys():
          if wet_dry_2dhist:
              json_fileout='hist2d_{}_{}_{:d}.json'.format(reg,
                                                              wet_dry_2dhist,
                                                              int(float(thrs)),)
          else:
              json_fileout='run_{}_{}.json'.format(reg, bintype)

          ### make sure the right python is called (python2.?)
          cmd = 'python2.7 process_data_europe.py -i {} -n {} -m {} -d {} -r masked -s {} -q {} -f {} -w {} -t {} -o {}'.format(
                                                                                runlist[res],
                                                                                res,
                                                                                maskedreglist[reg],
                                                                                bintype,
                                                                                ' '.join(se for se in season_list),
                                                                                s,
                                                                                frequency,
                                                                                wet_dry_2dhist,
                                                                                thrs,
                                                                                json_fileout)

          print cmd
          if not run_on_spice:
             arguments = para_args(cmd = cmd, msg = res)
             call_args.append(arguments)
          else:
             lotus_cmd = 'mqsub -W 06:00 {}'.format(cmd)  #mqsub -W before
             argument = para_args(cmd=lotus_cmd,
                                  msg='failed to submit to LOTUS')
             shellcmd(argument)

# if not run_on_spice:
#    print 'launching {} parallel process_data'.format(len(runlist)+1)
#    p = mp.Pool(len(runlist)+1)
#    p.map(shellcmd, call_args)

