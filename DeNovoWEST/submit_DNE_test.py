'''
@queenjobo @ksamocha 16/08/2019

Script that prints commands to run enrichment test - can adapt for computing environment
'''

__version__ = 1.1
__date__ = '2019-08-16'
__author__ = 'Joanna Kaplanis (@queenjobo) and Kaitlin Samocha (@ksamocha)'

import os
import argparse
import sys
import datetime


def get_rates_files(path):
    ratesfiles = [file for file in os.listdir(path) if (file.endswith(".txt.gz") and file.startswith('all_rates_'))]
    return(ratesfiles)


def submit_jobs(rates_path, weight_path, denovos_path, out_path, nmal, nfem, today_info):
    ratesfiles = get_rates_files(rates_path)
    for i,file in enumerate(ratesfiles):
        outfile = out_path + "dne_test_" + "_".join(file.split(".")[0].split("_")[2:4]) + "_ppv_{0}.tab".format(today_info)
        jobfile = out_path + str(i) + ".out"
        command = " 'python DNE_test.py --weightdic " + weight_path + ' --nmales ' + str(nmal) + ' --nfemales ' + str(nfem) + ' --denovos ' + denovos_path + ' --rates ' + rates_path + file + ' --output ' + outfile + "'"
        #uncomment out the following two lines if running on lsf
        #joboptions = "bsub -R'select[mem>20000] rusage[mem=20000]' -M20000 -o " + jobfile
        #command = joboptions + command
        print(command)

    
def main(args):
    current_date = datetime.date.today()
    today_info = '{0}_{1:02d}_{2:02d}'.format(current_date.year, current_date.month, current_date.day)
    submit_jobs(args.ratespath, args.weightspath, args.denovospath, args.outpath, args.nmale, args.nfemale, today_info)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
        Generate commands for DNE_test.py
    ''')
    parser.add_argument('--ratespath', action='store', type=str, dest='ratespath',
                        required=True, help='File path to rates files')
    parser.add_argument('--weightspath', action='store', type=str, dest='weightspath',
                        required=True, help='File path to where weights should be stored')
    parser.add_argument('--denovospath', action='store', type=str, dest='denovospath',
                        required=True, help='File path for the de novo variants')
    parser.add_argument('--outpath', action='store', type=str, dest='outpath',
                        required=True, help='Output file path')
    parser.add_argument('--nmale', action='store', type=int, dest='nmale',
                        required=True, help='Number of males')
    parser.add_argument('--nfemale', action='store', type=int, dest='nfemale',
                        required=True, help='Number of females')

    args = parser.parse_args()

    if not os.path.exists(args.ratespath):
        sys.exit('{0}: No such file or directory'.format(args.ratespath))
    if not os.path.exists(args.weightspath):
        sys.exit('{0}: No such file or directory'.format(args.weightspath))
    if not os.path.exists(args.denovospath):
        sys.exit('{0}: No such file or directory'.format(args.denovospath))

    main(args)
