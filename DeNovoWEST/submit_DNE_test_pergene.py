'''
@queenjobo @ksamocha 04/03/2019
Script that prints commands to run enrichment test

Adapted from submit_DNE_test.py to make temporary files per gene for the
de novo mutations in that gene and then to submit a job to the farm (or
print out the python commands).

Version 1.1: Modified to account for changes in DNE_test, specifically no
longer splitting the rates file. Also added a merge job at the end.
WARNING: Job array currently not working.

Version 2.0: Requires a missense weight dictionary and submits those jobs
as well.
Version 2.1: Instead of per gene files, prints 10 genes per file (by default).
Can modify the number of genes printed per file with --ngenesperfile.
Version 2.2: Fixing the --onlyprint option.
Version 2.3: Setting the queue to long and giving up on subprocess. Slight change
to the merge commands.
Version 2.4: Changing where the .out writes and the name of the output for the all
enrichment
Version 2.5: Adding --pvalcap flag to pass to the dne_test (float with default of 1)
Version 2.6: Avoiding hardcoding the header in the de novo file
Version 2.7: Fixing an error where an extra empty file would be printed when using
--ngenesperfile 1

'''

__version__ = 2.7
__date__ = '2020-03-23'
__author__ = 'Joanna Kaplanis (@queenjobo) and Kaitlin Samocha (@ksamocha)'

import os
import argparse
import sys
import datetime
import tempfile
import gzip


def submit_jobs(nfiles, temp_dir, rates_path, weight_path, mis_weight_path, out_path, nmal, nfem, today_info, pvalcap):
    '''Written to take a list of genes and submit DNE_test.py jobs based on their specific DNM and rates files

    Modified to use only one large rate file and the command:
    python /nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py --weightdic X --nmales X --nfemales X --denovos X --rates X --output X
    and potentially add --pvalcap X
    '''
    
    dnm_file = os.path.join(temp_dir, "tmp.\$LSB_JOBINDEX\.txt")

    outfile = os.path.join(temp_dir, "tmp.all.\$LSB_JOBINDEX\.output")
    misoutfile = os.path.join(temp_dir, "tmp.mis.\$LSB_JOBINDEX\.output")

    jobfile = os.path.join(temp_dir, 'pergene_all_dne_test_{0}.out'.format(today_info))
    misjobfile = os.path.join(temp_dir, 'pergene_mis_dne_test_{0}.out'.format(today_info))

    # make command
    all_command_array = ['\"', 'python', '/nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py',
                         '--weightdic', weight_path,
                         '--nmales', str(nmal),
                         '--nfemales', str(nfem),
                         '--denovos', dnm_file,
                         '--rates', rates_path,
                         '--output', outfile, '\"'] #,'--nsim', '1000', '\"']
    mis_command_array = ['\"', 'python', '/nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py',
                         '--weightdic', mis_weight_path,
                         '--nmales', str(nmal),
                         '--nfemales', str(nfem),
                         '--denovos', dnm_file,
                         '--rates', rates_path,
                         '--output', misoutfile, '\"'] # ,'--nsim', '1000', '\"'] 

    # add pvalcap
    if pvalcap != 1.0:
        all_command_array = all_command_array[:-1] + ['--pvalcap', str(pvalcap), '\"']
        mis_command_array = mis_command_array[:-1] + ['--pvalcap', str(pvalcap), '\"']
    
    # assuming LSF
    job_name = 'dne_test'
    job_id = '{0}[1-{1}]'.format(job_name, nfiles)
    joboptions_array = ['bsub', '-J', '\"{0}\"'.format(job_id), '-o', jobfile, '-q', 'long', '-M20000', '-R', '"select[mem>20000] rusage[mem=20000]"']

    # missense
    misjob_name = 'dne_test_mis'
    misjob_id = '{0}[1-{1}]'.format(misjob_name, nfiles)
    misjoboptions_array = ['bsub', '-J', '\"{0}\"'.format(misjob_id), '-o', misjobfile, '-q', 'long', '-M20000', '-R', '"select[mem>20000] rusage[mem=20000]"']
        
    full_command_array = joboptions_array + all_command_array
    fullmis_command_array = misjoboptions_array + mis_command_array
    print(' '.join(full_command_array))
    print(' '.join(fullmis_command_array))

    # submit merge job
    merge_id = 'dne_test_merge'
    merge_command = ["bash", "-c", "\"", "head", "-n", "1", os.path.join(temp_dir, "tmp.all.1.output"), ">",
                     '{0}merged_all_dne_test_ppv_{1}.tab'.format(out_path, today_info), "; tail", "-q", "-n", "+2",
                     os.path.join(temp_dir, "tmp.all.*.output"), "|", "sort", ">>",
                     '{0}merged_all_dne_test_ppv_{1}.tab'.format(out_path, today_info), "\""]
    mis_merge_command = ["bash", "-c", "\"", "head", "-n", "1", os.path.join(temp_dir, "tmp.mis.1.output"),
                         ">", '{0}merged_mis_dne_test_ppv_{1}.tab'.format(out_path, today_info), "; tail", "-q", "-n", "+2",
                         os.path.join(temp_dir, "tmp.mis.*.output"), "|", "sort", ">>",
                         '{0}merged_mis_dne_test_ppv_{1}.tab'.format(out_path, today_info), "\""]
    merge_joboptions = ['bsub', '-o', 'bsubout_merge_dne_test.out', '-J', '\"{0}\"'.format(merge_id),
                        '-M1000', '-R', '"select[mem>1000] rusage[mem=1000]"', '-w', '\"{0}\"'.format(job_id)]
        
    full_merge_job = merge_joboptions + merge_command
    fullmis_merge_job = merge_joboptions + mis_merge_command
    print(' '.join(full_merge_job))
    print(' '.join(fullmis_merge_job))


def print_jobs(nfiles, temp_dir, rates_path, weight_path, mis_weight_path, out_path, nmal, nfem, today_info, pvalcap):
    '''Prints the jobs to submit for a list of files

    Modified to use only one large rate file and the command:
    python /nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py --weightdic X --nmales X --nfemales X --denovos X --rates X --output X
    and maybe add --pvalcap X
    '''

    for i in range(1,nfiles+1):
        dnm_file = os.path.join(temp_dir, "tmp.{0}.txt".format(i))
        outfile = os.path.join(temp_dir, "tmp.all.{0}.output".format(i))
        misoutfile = os.path.join(temp_dir, "tmp.mis.{0}.output".format(i))

        all_command = 'python /nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py --weightdic {0} --nmales {1} --nfemales {2} --denovos {3} --rates {4} --output {5}'.format(weight_path, nmal, nfem, dnm_file, rates_path, outfile)
        mis_command = 'python /nfs/ddd0/jk18/dne/DeNovoWEST/DeNovoWEST/DNE_test.py --weightdic {0} --nmales {1} --nfemales {2} --denovos {3} --rates {4} --output {5}'.format(mis_weight_path, nmal, nfem, dnm_file, rates_path, misoutfile)

        # add pvalcap
        if pvalcap != 1.0:
            all_command = all_command + ' --pvalcap {0}'.format(pvalcap)
            mis_command = mis_command + ' --pvalcap {0}'.format(pvalcap)
        
        print(all_command)
        print(mis_command)

    # print the merge commands too
    print('### MERGE COMMANDS BELOW ###')
    merge_command = "bash -c \" head -n 1 {0} > {1} ; tail -q -n +2 {2} | sort >> {1} \" ".format(os.path.join(temp_dir, "tmp.all.1.output"), '{0}merged_all_dne_test_ppv_{1}.tab'.format(out_path, today_info), os.path.join(temp_dir, "tmp.all.*.output"))
    mis_merge_command = "bash -c \" head -n 1 {0} > {1} ; tail -q -n +2 {2} | sort >> {1} \" ".format(os.path.join(temp_dir, "tmp.mis.1.output"), '{0}merged_mis_dne_test_ppv_{1}.tab'.format(out_path, today_info), os.path.join(temp_dir, "tmp.mis.*.output"))

    print(merge_command)
    print(mis_merge_command)

    
def load_dnms_print(dnm_filepath, temp_dir, ngenesperfile):
    '''Print a temporary file for every gene in the de novo file separately

    Modified to print X genes per file, where X is defined by the user

    0: alt, 1: altprop_child, 2: chrom, 3: consequence, 4: cq, 5: hgnc_id, 6: id,
    7: maf, 8: pos, 9: prob, 10: raw, 11: ref, 12: score, 13: study, 14: symbol,
    15: constrained, 16: shethigh

    Now modified to deal with the fact that some de novo files have different headers
    0: symbol, 1: hgnc_id, 2: alt, 3: altprop_child, 4: chrom, 5: consequence, 6: cq, 7: id,
    8: maf, 9: pos, 10: prob, 11: raw, 12: ref, 13: score, 14: study, 15: constrained, 16: shethigh
    '''

    dnm_genes = {}
    my_open = gzip.open if (dnm_filepath.endswith('.gz') or dnm_filepath.endswith('.bgz')) else open

    with my_open(dnm_filepath, 'rt') as dnmfile:
        for line in dnmfile:
            line = line.split('\t')

            # header line -- save
            if 'symbol' in line:
                header = line
                # find symbol position
                symbol_col = line.index('symbol')
                continue

            # dnm line
            if line[symbol_col] in dnm_genes.keys():
                # add to entry
                dnm_genes[line[symbol_col]].append('\t'.join(line))
            else:
                dnm_genes[line[symbol_col]] = ['\t'.join(line)]

    sys.stderr.write('Writing files...\n')

    ngenes = len(dnm_genes.keys())
    counter_gene = 0
    counter_files = 1
    running_gene_counter = 1

    # open first file
    current_filehandle = temp_dir + '/tmp.{0}.txt'.format(counter_files)
    current_file = open(current_filehandle, 'w')
    current_file.write('\t'.join(header))
    
    # all dnms saved in dictionary now
    for gene_entry in dnm_genes.keys():
        # to keep track when running
        counter_gene += 1
        if (counter_gene % 1000) == 0:
            sys.stderr.write('\t{0} genes\n'.format(counter_gene))

        # If running counter < ngenesperfile
        if running_gene_counter < ngenesperfile:
            # write the dnms
            for dnm in dnm_genes[gene_entry]:
                current_file.write(dnm)
            running_gene_counter += 1
                
        # If running counter == ngenesperfile, write then close, open new
        elif running_gene_counter == ngenesperfile:
            # write the dnms
            for dnm in dnm_genes[gene_entry]:
                current_file.write(dnm)

            # close the current file
            current_file.close()

            # update the counters, reset running to 1
            counter_files += 1
            running_gene_counter = 1

            # if the counter_files is greater than the nuber of genes, break
            if counter_files > ngenes:
                counter_files -= 1
                break

            # open the new file and write the header
            current_filehandle = temp_dir + '/tmp.{0}.txt'.format(counter_files)
            current_file = open(current_filehandle, 'w')
            current_file.write('\t'.join(header))

        else:
            sys.exit('Somehow running_gene_counter is larger than the ngenesperfile limit\t running {0}\tngenes {1}'.format(running_gene_counter, ngenesperfile))

    # close the last file
    current_file.close()
    sys.stderr.write('Wrote {0} gene files\n'.format(counter_files))
    
    return(dnm_genes, counter_files)

    
def main(args):
    current_date = datetime.date.today()
    today_info = '{0}_{1:02d}_{2:02d}'.format(current_date.year, current_date.month, current_date.day)

    # Step 0: Make a temporary directory for DNM files
    # taken from Jeremy McRae's denovonear code (https://github.com/jeremymcrae/denovonear/blob/master/scripts/run_batch.py)
    temp_dir = tempfile.mkdtemp(dir=args.outpath)

    # Step 1: load all DNMs and print temporary files per gene
    sys.stderr.write('Loading DNMs...\n')
    (dnm_genes, nfiles) = load_dnms_print(args.denovospath, temp_dir, args.ngenesperfile)
    sys.stderr.write('There are {0} genes with DNMs\n'.format(len(dnm_genes.keys())))

    # Step 2: Submit jobs
    sys.stderr.write('Preparing jobs to submit...\n')
    if args.onlyprint:
        # separate print script
        print_jobs(nfiles, temp_dir, args.ratespath, args.weightspath, args.misweightspath,
                   args.outpath, args.nmale, args.nfemale, today_info, args.pvalcap)
    else:
        submit_jobs(nfiles, temp_dir, args.ratespath, args.weightspath, args.misweightspath,
                    args.outpath, args.nmale, args.nfemale, today_info, args.pvalcap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
        Generate commands per gene for DNE_test.py
    ''')
    parser.add_argument('--ratespath', action='store', type=str, dest='ratespath',
                        required=True, help='File path to rates file')
    parser.add_argument('--weightspath', action='store', type=str, dest='weightspath',
                        required=True, help='File path to weights file')
    parser.add_argument('--misweightspath', action='store', type=str, dest='misweightspath',
                        required=True, help='File path to missense weights file')
    parser.add_argument('--denovospath', action='store', type=str, dest='denovospath',
                        required=True, help='File path for the de novo variants')
    parser.add_argument('--outpath', action='store', type=str, dest='outpath',
                        required=True, help='Output file path')
    parser.add_argument('--nmale', action='store', type=int, dest='nmale',
                        required=True, help='Number of males')
    parser.add_argument('--nfemale', action='store', type=int, dest='nfemale',
                        required=True, help='Number of females')

    # optional flag to only print the command -- not submit
    parser.add_argument('--onlyprint', action='store_true', dest='onlyprint',
                        help='Only print the commands, do not submit')
    # optional flag to change the default number of genes per file from 10 to X
    parser.add_argument('--ngenesperfile', action='store', type=int, dest='ngenesperfile',
                        default = 10, help='Number of genes per output file')
    # optional flag to change the p-value cap in DNE test
    parser.add_argument('--pvalcap', action='store', type=float, dest='pvalcap',
                        default=1.0, help='Cap for the p-value in DNE test (default is 1.0)')

    args = parser.parse_args()

    if not os.path.exists(args.ratespath):
        sys.exit('{0}: No such file or directory'.format(args.ratespath))
    if not os.path.exists(args.weightspath):
        sys.exit('{0}: No such file or directory'.format(args.weightspath))
    if not os.path.exists(args.misweightspath):
        sys.exit('{0}: No such file or directory'.format(args.misweightspath))
    if not os.path.exists(args.denovospath):
        sys.exit('{0}: No such file or directory'.format(args.denovospath))

    main(args)
