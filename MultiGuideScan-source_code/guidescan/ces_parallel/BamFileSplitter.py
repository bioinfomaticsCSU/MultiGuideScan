__author__ = 'Alex Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Split BAM files by a set amount of lines

bam_file_split(bam_dir,bam_name,split_lines,prefix_name) = splits a BAM file into subset files containing at most
														   a 'split_lines' amount of lines. Each subset file contains
														   the prefix which is defined in 'prefix_name'. Subset files
														   are located in 'bam_dir'. The exact BAM file that is used is
														   defined by bam_dir/bam_name.
"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import argparse

#############################
#                           #
#   Auxillary Functions     #
#                           #
#############################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--dir','-d',help='filepath to directory which hosts BAM file',required=True)
	parser.add_argument('--name','-n',help='name of BAM file',required=True)
	parser.add_argument('--lines','-l',help='amount of lines in each split file',required=True)
	parser.add_argument('--prefix','-p',help='prefix added to the start of each split filename',required=True)

	args = parser.parse_args()
	bam_dir = args.dir
	bam_name = args.name
	split_lines = args.lines
	prefix_name = args.prefix

	return bam_dir,bam_name,split_lines,prefix_name

def bam_file_split(bam_dir,bam_name,split_lines,prefix_name):
	"""Split a bam file by lines

	Convention: A bam file is the input to this function and is split into smaller bam files by the number of lines
	present. The function takes as input the filepath to the directory which hosts the bam file, the bam filename, the
	amount of lines from the original bam file the user wants each split file to contain, and the prefix associated with
	each split file. The split files are written to the bam_dir directory, consequently segment_name should be distinct
	from the filename of the original file. The function relies on system calls to Unix and Samtools. Samtools should
	be in the user's PATH. The function writes out an intermediate SAM file entitled header.sam which is deleted at the
	end of the function's run.

	Args:
	bam_dir: the filepath to the directory which hosts the bam file to be split

	bam_name: the filename of the bam file to be split

	split_lines: the amount of lines each split file will contain from the original bam file

	prefix_name: prefix attached to the start of each split file.
	"""

	#get bam filepath
	bam_file = bam_dir + '/' + bam_name

	#split bam file
	os.chdir(bam_dir)
	cmd1 = 'samtools view ' + bam_file +' | split -l ' + str(split_lines) + ' ' + '-' + ' ' + prefix_name
	os.system(cmd1)

	#get all the split files
	split_files = []
	for File in os.listdir(bam_dir):
		if File.startswith(prefix_name) and not File.endswith('.bam'):
			split_files.append(File)

	#get header
	cmd2 = 'samtools view -H ' + bam_file + ' > header.sam'
	os.system(cmd2)

	#recreate bam files with header
	for i in range(len(split_files)):
		cmd3 = 'cat header.sam ' + split_files[i] + ' | samtools view -Sb - > ' + split_files[i] + '.bam'
		cmd4 = 'rm ' + split_files[i]
		cmd5 = 'samtools index ' + split_files[i] + '.bam'
		os.system(cmd3)
		os.system(cmd4)
		os.system(cmd5)

	cmd5 = 'rm header.sam'
	os.system(cmd5)

	return 'bam file splitting complete'

#####################
#                   #
#   Main Function   #
#                   #
#####################

if __name__ == '__main__':
	bam_dir,bam_name,split_lines,prefix_name = arg_parser()
	msg1 = bam_file_split(bam_dir,bam_name,split_lines,prefix_name)
	print msg1
	print 'Finished'




