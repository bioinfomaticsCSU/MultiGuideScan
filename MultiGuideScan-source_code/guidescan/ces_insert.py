__author__ = 'Alex Perez and Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""This script will take cutting efficiency scores and insert them into a BAM file as a tag
"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import subprocess
import sys
import argparse

#############################
#                           #
#   Auxillary Functions     #
#                           #
#############################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--dir','-d',help='filepath to the directory hosting the score file and database BAM file',
						required=True)
	parser.add_argument('--bam','-b',help='name of database BAM file',required=True)
	parser.add_argument('--score','-s',help='name of file with cutting efficiency scores',required=True)

	args = parser.parse_args()
	file_dir = args.dir
	bam_filename = args.bam
	score_filename = args.score

	return file_dir,bam_filename,score_filename

def header_writeout_and_strip(file_dir,filename):
	"""extracts header from BAM file and converts BAM to SAM

	Args:
	file_dir = filepath to directory hosting the BAM file
	filename = name of BAM file containing sgRNAs
	"""
	cmd1 = 'samtools view -H ' + str(file_dir) + '/' + str(filename) + ' > ' + str(file_dir) + '/' + 'header.txt'
	cmd2 = "samtools view " + file_dir + '/' + filename + " | awk '!/@/' - " + ' > ' + file_dir + '/' + 'noheader_' + filename
	os.system(cmd1)
	os.system(cmd2)

	return ('header written out and stripped from %s \n' % (filename))

def files_of_equal_length(file_dir,bam_filename,score_filename):
	"""check score file and SAM file have same amount of lines

	convention: this is a sanity check to ensure that before any heavy computation occurs
	that the amount of lines in the scores file and the amount of lines in the sgRNA BAM
	file are the same. There needs to be a one-to-one concordance since scores are going to
	be added directly to the BAM file. If there is not one-to-one concordance than a system
	exit occurs. The bam_filename will be contained within the score_filename. This will occur
	automatically when scores are computed using CutEffScore.py.

	Args:
	file_dir = filepath to directory hosting the BAM file and scores file
	bam_filename = name of BAM file containing sgRNAs
	score_filename = name of scores file containg sgRNA scores.
	"""
	cmd1 = 'wc -l ' + file_dir + '/noheader_' + bam_filename
	cmd2 = 'wc -l ' + file_dir + '/' + score_filename

	out1 = subprocess.check_output(cmd1,shell=True)
	out2 = subprocess.check_output(cmd2,shell=True)
	out1_val = int(out1.lstrip().rstrip().split()[0])
	out2_val = int(out2.lstrip().rstrip().split()[0])

	if out1_val == out2_val:
		return ('%s and %s are of equal length \n' % (bam_filename,score_filename))
	else:
		sys.stderr.write('%s and %s are of unequal length \n' % (bam_filename,score_filename))
		sys.exit()

def BAM_file_modify(file_dir,bam_filename,score_filename):
	"""add cutting efficiency scores to BAM file

	Convention: this function make a series of system calls using UNIX and AWK commands. The
	function takes the sgRNA BAM file and scores file and first checks to see if they are sorted
	in the same order. It does this by looking at the order of the unique identifiers in both files
	(sgRNA target sequence). If the files are in the same order then cutting efficiency scores are
	added to the BAM file. If the files are not sorted in the same order then both files are sorted
	by their unique identifiers. Before scores are added to the BAM file a check is done to make sure
	that the elements in the two sorted files are identical. If the elements and their order are not
	identicial a system exit occurs, otherwise the scores are added to the BAM file.

	Args:
	file_dir = filepath to directory hosting the BAM file and scores file
	bam_filename = name of BAM file containing sgRNAs
	score_filename = name of scores file containg sgRNA scores.
	"""
	cmd1 = "awk '{print $1}' " + file_dir + '/' + 'noheader_' + bam_filename + ' > ' + file_dir + '/bam_sequences_' + bam_filename
	cmd2 = "awk '{print $1}' " + file_dir + '/' + score_filename + ' > ' + file_dir + '/score_sequences_' + score_filename
	cmd3 = 'diff ' + file_dir + '/bam_sequences_' + bam_filename + ' ' + file_dir + '/score_sequences_' + score_filename + ' > ' + file_dir + '/bam_score_diff_' + bam_filename
	cmd4 = 'wc -l ' + file_dir + '/bam_score_diff_' + bam_filename
	cmd5 = 'rm ' + file_dir + '/bam_sequences_' + bam_filename
	cmd6 = 'rm ' + file_dir + '/score_sequences_' + score_filename
	cmd7 = 'rm ' + file_dir + '/bam_score_diff_' + bam_filename
	cmd8 = "awk '{print $6}' " + file_dir + '/' + score_filename + ' > ' + file_dir + '/scores_' + score_filename
	cmd9 = 'paste ' + file_dir + '/' + 'noheader_' + bam_filename + ' ' + file_dir + '/scores_'+score_filename+ ' > ' + file_dir + '/CutEffScore_' + bam_filename
	cmd10 = 'rm ' + file_dir + '/scores_' + score_filename
	cmd11 = 'sort -k1,1 ' + file_dir + '/' + 'noheader_' + bam_filename + ' > ' + file_dir + '/' + 'noheader_A' + bam_filename
	cmd12 = 'sort -k1,1 ' + file_dir + '/' + score_filename + ' > ' + file_dir + '/A' + score_filename
	cmd13 = "awk '{print $6}' " + file_dir + '/A' + score_filename + ' > ' + file_dir + '/scores_' + score_filename
	cmd14 = 'paste ' + file_dir + '/' + 'noheader_A' + bam_filename + ' ' + file_dir + '/scores_'+score_filename+' > ' + file_dir + '/CutEffScore_' + bam_filename
	cmd15 = 'rm ' + file_dir + '/' + 'noheader_A' + bam_filename
	cmd16 = 'rm ' + file_dir + '/A' + score_filename
	cmd17 = "awk '{print $1}' " + file_dir + '/' + 'noheader_A' + bam_filename + ' > ' +  file_dir + '/bamSortConfirm' + bam_filename
	cmd18 = "awk '{print $1}' " + file_dir + '/A' + score_filename + ' > ' +  file_dir + '/ScoreSortConfirm' + score_filename
	cmd19 = 'diff ' + file_dir + '/bamSortConfirm' +bam_filename+ ' ' + file_dir + '/ScoreSortConfirm'+score_filename + ' > ' + file_dir + '/SortConfirm_diff'+bam_filename
	cmd20 = 'wc -l ' + file_dir + '/SortConfirm_diff'+bam_filename
	cmd21 = 'rm ' + file_dir + '/bamSortConfirm'+bam_filename
	cmd22 = 'rm ' + file_dir + '/ScoreSortConfirm'+score_filename
	cmd23 = 'rm ' + file_dir + '/SortConfirm_diff'+bam_filename

	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	result = subprocess.check_output(cmd4,shell=True)
	result_val = int(result.lstrip().rstrip().split()[0])
	os.system(cmd5)
	os.system(cmd6)
	os.system(cmd7)
	if result_val == 0:
		sys.stdout.write('%s and %s lines are in identical order \n' % (bam_filename,score_filename))
		os.system(cmd8)
		os.system(cmd9)
		os.system(cmd10)
		sys.stdout.write('%s adjusted to contain cutting efficiency scores \n' % (bam_filename))
	else:
		sys.stdout.write('%s and %s lines are not in identical order: files will now be sorted \n' %
						 (bam_filename,score_filename))
		os.system(cmd11)
		os.system(cmd12)
		os.system(cmd17)
		os.system(cmd18)
		os.system(cmd19)
		result2 = subprocess.check_output(cmd20,shell=True)
		result_val2 = int(result2.lstrip().rstrip().split()[0])
		os.system(cmd21)
		os.system(cmd22)
		os.system(cmd23)
		if result_val2 == 0:
			os.system(cmd13)
			os.system(cmd14)
			os.system(cmd10)
			os.system(cmd15)
			os.system(cmd16)
			sys.stdout.write('%s sorted and adjusted to contain cutting efficiency scores \n' % (bam_filename))
		else:
			sys.stderr.write('ERROR: sorted files, %s and %s do not contain identical elements \n' %
							 (bam_filename,score_filename))
			sys.exit()
	return ('%s modified to include cutting efficiency tag \n' % (bam_filename))

def regenerate_BAM_file(file_dir,bam_filename):
	"""converts SAM file with cutting efficiency scores back to BAM file

	Args:
	file_dir = filepath to directory hosting the BAM file and scores file
	bam_filename = name of BAM file containing sgRNAs

	"""
	cmd1 = 'cat ' + file_dir + '/header.txt ' + file_dir + '/CutEffScore_' + bam_filename + ' > ' + file_dir + '/' + 'versionA_' + bam_filename
	cmd2 = 'samtools view -Sb ' + file_dir + '/versionA_' + bam_filename + ' > ' + file_dir + '/modified_' + bam_filename
	cmd3 = 'samtools sort ' + file_dir + '/modified_' + bam_filename + ' -o ' + file_dir + '/cutting_efficiency_scores_added_' + bam_filename
	cmd4 = 'rm ' + file_dir + '/versionA_' + bam_filename
	cmd5 = 'rm ' + file_dir + '/modified_' + bam_filename
	cmd6 = 'rm ' + file_dir + '/CutEffScore_' + bam_filename

	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	os.system(cmd4)
	os.system(cmd5)
	os.system(cmd6)
	return 'database bam file regenerated \n'

def clean_up(file_dir,bam_filename):
	"""clean up intermediate files
	"""
	cmd1 = 'rm header.txt'
	cmd2 = 'rm ' + file_dir + '/noheader_' + bam_filename

	#os.system(cmd1)
	os.system(cmd2)
	return 'clean up complete \n'

def BAM_index(file_dir,bam_filename):
	cmd = 'samtools index ' + file_dir + '/cutting_efficiency_scores_added_' + bam_filename
	os.system(cmd)

def file_path_verification(database):
	"""function that verifies that the provided filepath exists before further processing occurs

	Input:
	database: absolute filepath to the sgRNA database
	"""
	filename = database.split('/')[-1]
	if os.path.exists(database):
		sys.stdout.write('%s located \n' % (filename))
		return 0
	else:
		sys.stderr.write('ERROR: %s does not exist at the filepath indicated: please provide correct absolute filepath \n'
						 % (filename))
		return 1
		#sys.exit()

#####################
#                   #
#   Core Function   #
#                   #
#####################

def core(file_dir,bam_filename,score_filename):

	#ensure BAM file and score files exist before computation
	bam_exists = file_path_verification(file_dir + '/' + bam_filename)
	score_exists = file_path_verification(file_dir + '/' + score_filename)

	if bam_exists or score_exists != 0:
		sys.exit()

	else:
		os.chdir(file_dir)

		#convert BAM file to SAM file, keep header
		msg = header_writeout_and_strip(file_dir,bam_filename)
		sys.stdout.write(msg)

		#check to make sure scores file and SAM file are same length
		msg = files_of_equal_length(file_dir,bam_filename,score_filename)
		sys.stdout.write(msg)

		#add cutting efficiency scores to SAM file
		msg = BAM_file_modify(file_dir,bam_filename,score_filename)
		sys.stdout.write(msg)

		#recreate BAM file
		msg = regenerate_BAM_file(file_dir,bam_filename)
		sys.stdout.write(msg)

		#clean up
		msg = clean_up(file_dir,bam_filename)
		sys.stdout.write(msg)

		#index BAM file
		BAM_index(file_dir,bam_filename)

		#print 'Finished'
		sys.stdout.write('cutting efficiency scores inserted into %s \n' % (bam_filename))
		return
