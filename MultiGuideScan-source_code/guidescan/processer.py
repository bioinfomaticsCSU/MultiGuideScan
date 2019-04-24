__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Master script for guideRNA database construction.

Save user input parameters and run all subroutines.

Subroutines are:
extract k-mers (candidate guideRNAs) from FASTA
count k-mer occurrences
build a trie with all k-mers
label k-mers that cannot be guideRNAs
collect info about off-targets of all remaining candidate guideRNAs
    and produce guideRNA database in BAM format

Then access guideRNA database using guidequery.py.

"""

#################
#               #
#   Libraries   #
#               #
#################

import argparse
import os

import util
import kmers
import guides
import bamdata

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
    p = argparse.ArgumentParser(description='Precompute guideRNA database from'
                                            ' genomic sequences',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('-f', dest='fasta', required=True,
                   help='path to fasta file or folder with fasta files'
                        ' (will use all .fa, .fasta, .fa.gz, .fasta.gz files'
                        ' found in the folder)')
    p.add_argument('-n', dest='name', default='myguides',
                   help='project name, use in all output (will produce'
                        ' a folder with this name containing intermediate and'
                        ' final files in it)')
    p.add_argument('--minchr', dest='minchrlen', type=int, default=10000,
                   help='minimum chromosome length to consider,'
                        ' chromosomes in input FASTA that are shorter than'
                        ' this will be excluded from analysis;'
                        ' simple way to exclude scaffolds unassigned to known'
                        ' chromosomes etc.')
    p.add_argument('-c', dest='chrom', default='',
                   help='list names of chromosomes (comma separated) that will'
                        ' be used in analysis, or name of file where this list'
                        ' is stored')
    p.add_argument('-l', dest='length', type=int, default=20,
                   help='desired length of guideRNAs (not including PAM)')
    p.add_argument('-p', dest='pam', default='NGG',
                   help='PAM sequence')
    p.add_argument('-a', dest='altpam', default='NAG',
                   help='alternative PAM sequences (separate multiple ones'
                        ' by commas), will not be used in primary guideRNAs,'
                        ' but will be considered in off-targets;'
                        ' all PAM sequences should be mutually exclusive'
                        ' and of the same length')
    p.add_argument('--pampos', dest='pampos', default='end',
                   choices=['start', 'end'],
                   help='position of PAM with respect to guideRNA')
    p.add_argument('-s', dest='sim', type=int, default=2,
                   help='minimum mismatch similarity between guideRNAs;'
                        ' a candidate guideRNA (with primary PAM)'
                        ' should not have alternative occurences in the genome'
                        ' (with any PAM) with less than'
                        ' this many mismatches (not including PAM)')
    p.add_argument('-d', dest='offdist', type=int, default=3,
                   help='maximum distance to consider from guideRNA to its'
                        ' off-target;'
                        ' off-target is an alternative occurrence (with any'
                        ' PAM) of this guideRNA in the genome at edit distance'
                        ' at most this number (including PAM);'
                        ' currently values larger than 4 are infeasible for'
                        ' large (e.g., mammalian) genomes, and value 3'
                        ' will take long time to compute; use -1 if do not'
                        ' want any off-target info in resulting database'
                        ' (can add it later using bamdata)')
    p.add_argument('-k', dest='greateroffdist', type=int, default=4,
                   help='a number greater than offdist used for preprocessed data'
                        '(the length of key for classifying guide RNAs)')
    p.add_argument('--maxoffpos', dest='maxoffpos', type=int, default=10,
                   help='maximum number of positions of k-mers to remember;'
                        ' for k-mer occurring multiple times in the genome '
                        ' (such k-mers cannot be guideRNAs, but their'
                        ' positions can be off-targets of guideRNAs)'
                        ' store at most this many arbitrary their occurrences'
                        ' in the genome')
    p.add_argument('--maxoffcount', dest='maxoffcount', type=int, default=1000,
                   help='maximum number of off-targets to store for'
                        ' a guideRNA in a resulting BAM library')
    p.add_argument('--gnupath', dest='gnupath', default='',
                   help='path to gnu utilities, e.g. "/usr/local/bin";'
                        ' if empty, use system defaults;'
                        ' requires: cut, sort, uniq, shuf')
    p.add_argument('-t', dest='threads', type=int, default=1,
                   help='how many threads to use; do not specify more'
                        ' than you have on your system;'
                        ' currently not implemented')

    args = p.parse_args()

    return args

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

    #user inputs
    args = arg_parser()
    args_dict = args.__dict__
    #tidy PAM and chrom args
    args_dict['altpam'] = [s.upper() for s in args_dict['altpam'].split(',')]
    args_dict['altpam'] = [s.strip() for s in args_dict['altpam'] if s]
    args_dict['pam'] = args_dict['pam'].upper()
    if args_dict['chrom']:
        if os.path.isfile(args_dict['chrom']):
            chroms = open(args_dict['chrom']).read().split(',')
        else:
            chroms = args_dict['chrom'].split(',')
        chroms = [c.strip() for c in chroms]
        chroms = [c for c in chroms if c]
    else:
        chroms = []
    args_dict['chrom'] = chroms

    util.print_log('save arguments...')
    util.print_args(args_dict)
    util.save_args(args_dict)
    util.print_log('done')

    #main
    util.print_log2('start extract_process_kmers()')
    kmers.extract_process_kmers(args_dict['name'])

    util.print_log2('start analyze_guides()')
    kmers_trie = guides.analyze_guides(args_dict['name'])

    util.print_log2('start produce_bams_main()')
    bamdata.produce_bams_main(kmers_trie, args_dict['name'])

    util.print_log2('processer done.')

if __name__ == '__main__':
    main()
