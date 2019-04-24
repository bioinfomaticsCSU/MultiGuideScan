__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""this script will insert cutting efficiency scores into sgRNA databases if scores and database are in same order
"""

#################
#               #
#   Libraries   #
#               #
#################

import argparse
import os
import sys
import subprocess

#############################
#                           #
#   Auxillary Functions     #
#                           #
#############################

def arg_parse():
	parser = argparse.ArgumentParser()
	parser.add_argument('--out','-o',help='output directory',required=True)
	parser.add_argument('--newDatabase','-n',help='full filepath to the new database (does not have scores)',required=True)
	parser.add_argument('--oldDatabase','-v',help='full filepath to the old database (has scores)',required=True)

	args = parser.parse_args()
	outdir = args.out
	db1 = args.newDatabase
	db2 = args.oldDatabase

	return outdir,db1,db2

def BAM_index(file_dir,bam_filename):
	cmd = 'samtools index ' + file_dir + '/' + bam_filename
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

def core(indir1,db1,indir2,db2):

	#user inputs
	#outdir,db1,db2 = arg_parse()

	#ensure paths exist
	indir1_exist = file_path_verification(indir1)
	indir2_exist = file_path_verification(indir2)
	db1_exist = file_path_verification(indir1 + '/' + db1)
	db2_exist = file_path_verification(indir2 + '/' + db2)

	if indir1_exist or indir2_exist or db1_exist or db2_exist != 0:
		sys.exit(1)

	else:
		#generate database filepaths
		db1 = indir1 + '/' + db1
		db2 = indir2 + '/' + db2

		#ensure that scores and bam file have the same content and are in the same order
		cmd1 = "samtools view " + db1 + " | awk '{print $1}' - > " +indir1+"/testA.txt"
		cmd2 = "samtools view " + db2 + " | awk '{print $1}' - > " +indir1+"/testB.txt"
		cmd3 = 'diff testA.txt testB.txt > testC.txt'
		cmd4 = 'wc -l testC.txt'
		cmd5 = 'rm testA.txt'
		cmd6 = 'rm testB.txt'
		cmd7 = 'rm testC.txt'

		#get cutting efficiency scores from db2 file
		cmd8 = "samtools view " + db2 + " | awk '{print $NF}' - > " +indir1+"/cutting_efficinecy_scores.txt"

		#add cutting efficiency scores to new database
		cmd9 = "samtools view " + db1 + " | paste - cutting_efficinecy_scores.txt > newdatabase_w_scores.sam"
		cmd10 = 'rm cutting_efficinecy_scores.txt'
		cmd11 = "samtools view -H " + db1 + " > header.txt"
		cmd12 = "cat header.txt newdatabase_w_scores.sam > newdatabase_w_scores_wheader.sam"
		cmd13 = 'rm newdatabase_w_scores.sam'
		cmd14 = 'rm header.txt'
		cmd15 = "samtools view -bS newdatabase_w_scores_wheader.sam > newdatabase_w_transfered_scores.bam"
		cmd16 = 'rm newdatabase_w_scores_wheader.sam'

		#operating system calls
		os.chdir(indir1)
		os.system(cmd1)
		os.system(cmd2)
		os.system(cmd3)
		result = subprocess.check_output(cmd4,shell=True)
		result_val = int(result.lstrip().rstrip().split()[0])
		os.system(cmd5)
		os.system(cmd6)
		os.system(cmd7)
		if result_val == 0:
			os.system(cmd8)
			os.system(cmd9)
			os.system(cmd10)
			os.system(cmd11)
			os.system(cmd12)
			os.system(cmd13)
			os.system(cmd14)
			os.system(cmd15)
			os.system(cmd16)
		else:
			sys.stderr.write('BAM file and scoring file either do not contain same elements or not in same order \n')
			sys.exit(1)


		sys.stdout.write('cutting efficiency score transfer from %s to %s complete \n' % (db2,db1))

		BAM_index(indir1,'newdatabase_w_transfered_scores.bam')
		return

