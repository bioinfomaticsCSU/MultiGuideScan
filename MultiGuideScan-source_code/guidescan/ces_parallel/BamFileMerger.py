__author__ = 'Alex Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Merge split bam files into consolidated bam file

bam_file_merge(bam_dir,prefix_name,merged_file_name) = takes split bam files which are located in 'bam_dir' and each
                                                       has the prefix defined in 'prefix_name'. The resulting output
                                                       file is located in 'bam_dir' and is entitled 'merged_file_name'
                                                       it is the merged output of the split files.
"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import sys
import argparse

#############################
#                           #
#   Auxillary Functions     #
#                           #
#############################

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir','-d',help='filepath to the directory hosting the split bam files',required=True)
    parser.add_argument('--prefix','-p',help='name of the prefix of the split bam files',required=True)
    parser.add_argument('--merged','-m',help='name of the merged file',required=True)
    parser.add_argument('--score','-s',help='type of score, either specificity (specificity) or efficiency (efficiency)',required=True)

    args = parser.parse_args()
    bam_dir = args.dir
    prefix_name = args.prefix
    merged_file_name = args.merged
    score = args.score

    return bam_dir,prefix_name,merged_file_name,score

def bam_file_merge(bam_dir,prefix_name,merged_file_name,score):
    """Merge split bam files into consolidated bam file

    Convention: A series of split bam files are taken as input. The series of bam file are found in 'bam_dir'. The bam
    files which will be merged all have the prefix defined in 'prefix_name'. The resulting merged file will be found in
    'bam_dir' and will be entitled 'merged_file_name'.

    Args:
    bam_dir: the directory which hosts the split bam files to be merged

    prefix_name: the prefix associated with the split bam files

    merged_file_name: the name of the merged bam file resulting from the merging of split bam files
    """
    if score == 'specificity':
        file_string_prefix = 'cutting_specificity_scores_added_'
    elif score == 'efficiency':
        file_string_prefix = 'cutting_efficiency_scores_added_'
    else:
        sys.stderr.write('ERROR: score value %s is not accepted; enter either "specificity" or "efficiency" \n')
        sys.exit(1)

    os.chdir(bam_dir)
    cmd1 = 'samtools merge ' + merged_file_name
    for File in os.listdir(bam_dir):
        if File.startswith(file_string_prefix + prefix_name) and File.endswith('.bam'):
            cmd1 = cmd1 + ' ' + File

    #print cmd1
    os.system(cmd1)
    return 'bam files merged'

#####################
#                   #
#   Main Function   #
#                   #
#####################

if __name__ == '__main__':

    bam_dir,prefix_name,merged_file_name,score = arg_parser()
    msg1 = bam_file_merge(bam_dir,prefix_name,merged_file_name,score)
    print msg1
    print 'Finished'






