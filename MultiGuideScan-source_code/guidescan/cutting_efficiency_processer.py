__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""computes and inserts Rule Set 2 cutting efficiency scores into sgRNA database"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import os
import argparse

from pkg_resources import resource_exists, resource_filename

import ces_compute
import ces_insert
import ces_database_to_database_transfer

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d','--database_directory',help='absolute filepath to directory hosting sgRNA database with '
	                                                     'no cutting efficiency scores',required=True)
	parser.add_argument('-n','--database_name',help='filename of sgRNA database with no cutting efficiency scores',
	                    required=True)
	parser.add_argument('-f','--fasta_filepath',help='absolute filepath to organism FASTA file. Index for FASTA should '
													 'also be present in same directory (.fai). Must be single aggregate FASTA file for organism')
	parser.add_argument('-a','--amino',default=-1,help='amino acid information to be used in cutting efficiency.'
													   ' Default -1 indicating to not use amino acid information')
	parser.add_argument('-p','--peptide',default=-1,help='peptide information to be used in cutting efficiency.'
														 ' Default -1 indicating to not use peptide information')
	parser.add_argument('-d2','--database_directory2',help='absolute filepath to directory hosting sgRNA database with cutting '
														   'efficency scores already included and generated with all the '
														   'same run parameters as another database except for -d. This '
														   'parameter, coupled with -n2, allows for the transfer of '
														   'cutting efficiency scores between two sgRNA databases that '
														   'differ only in how many mismatches off-targets are '
														   'enumerated to')
	parser.add_argument('-n2','--database_name2',help='filename of sgRNA database with cutting efficiency scores already'
													  ' included and generated with all the same run parameters as '
													  'another database except for -d. This parameter, coupled with -d2,'
													  ' allows for the transfer of cutting efficiency scores between two'
													  ' sgRNA databases that differ only in how many mismatches '
													  'off-targets are enumerated to')

	args = parser.parse_args()
	indir = args.database_directory
	filename = args.database_name
	fasta = args.fasta_filepath
	aa_cut = args.amino
	per_peptide = args.peptide
	indir2 = args.database_directory2
	filename2 = args.database_name2

	return indir,filename,fasta,aa_cut,per_peptide,indir2,filename2

def database_check(database):
	"""function that verifies that a sgRNA database meets the requirements for Rule Set 2 on-target cutting efficiency
		score computation

		Input:
		database: output of GuideScan.processor. This will be the BAM file database of sgRNAs that meet the user's
		predefined parameters. The check is done by ensuring proper run parameters for Rule Set 2 scoring are present in
		bam header.

		Note:
		For a database to be eligible for Rule Set 2 on-target cutting efficiency scoring, it must be composed of sgRNAs
		whose complementary length region is 20bp, has only the NGG PAM as the canonical PAM (altpam can have any set of
		PAM sequences though), and the NGG PAM must be at the end (3' end) of the sgRNA target region.
	"""
	database_name = database.split('/')[-1]
	cmd = 'samtools view -H ' + database + ' | grep @CO > header_' + database_name.replace('.bam', '.txt')
	cmd2 = 'rm header_' + database_name.replace('.bam', '.txt')
	os.system(cmd)

	with open('header_' + database_name.replace('.bam', '.txt')) as f:
		exit_code = 0
		empty_header_pampos,empty_header_length,empty_header_pam = 0,0,0
		for line in f:
			clean_line = line.lstrip().rstrip()
			parts = clean_line.split()
			for i in range(len(parts)):
				if parts[i] == "'pampos':":
					if "'end'," == parts[i + 1] and len(parts[i + 1]) == 6:  # 6 because value is string of 6 characters: 'end',
						empty_header_pampos += 1
					else:
						exit_code += 1
						sys.stderr.write('ERROR: PAM is not at end (3 prime) end of sgRNA in %s \n' % (database_name))
				if parts[i] == "'length':":
					if "20," == parts[i + 1] and len(parts[i + 1]) == 3:  # 3 because value is string of 3 characters: 20,
						empty_header_length += 1
					else:
						exit_code += 1
						sys.stderr.write('ERROR: sgRNA complementary sequence length not 20bp in %s \n' % (database_name))
				if parts[i] == "'pam':":
					if "'NGG'," == parts[i + 1] and len(parts[i + 1]) == 6:  # 6 because value is string of 6 characters: 'NGG',
						empty_header_pam += 1
					else:
						exit_code += 1
						sys.stderr.write('ERROR: PAM sequence is not uniquely NGG in %s \n' % (database_name))
		os.system(cmd2)
		if exit_code != 0:
			sys.stderr.write('ERROR: database %s does not meet requirements for Rule Set 2 cutting efficiency scoring \n'
							 % (database_name))
			sys.exit(1)
		elif empty_header_pampos == 0:
			sys.stderr.write('ERROR: database %s does not have pampos field in header, computation cannot proceed \n'
							 % (database_name))
			sys.exit(1)
		elif empty_header_length == 0:
			sys.stderr.write('ERROR: database %s does not have length field in header, computation cannot proceed \n'
							 % (database_name))
			sys.exit(1)
		elif empty_header_pam == 0:
			sys.stderr.write('ERROR: database %s does not have pam field in header, computation cannot proceed \n'
							 % (database_name))
			sys.exit(1)
		else:
			sys.stdout.write('MSG: database %s meets requirements for Rule Set 2 cutting efficiency score computation \n'
							 % (database_name))
			return

def database_path_verification(database):
	"""function that verifies that the provided filepath exists before further processing occurs

	Input:
	database: absolute filepath to the sgRNA database
	"""
	filename = database.split('/')[-1]
	if os.path.exists(database):
		sys.stdout.write('%s located \n' % (filename))
		return
	else:
		sys.stderr.write('ERROR: %s does not exist at the filepath indicated: please provide correct absolute filepath \n'
						 % (filename))
		sys.exit(1)

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

	indir,filename,fasta,aa_cut,per_peptide,indir2,filename2 = arg_parser()

	#ensure sklearn is installed and required version
	try:
		import sklearn
		if sklearn.__version__ == '0.16.1':
			sys.stdout.write('sklearn version 0.16.1 identified \n')
			pass
		else:
			sys.stderr.write('cutting efficiency score model requires sklearn version 0.16.1, please install from %s \n' %
							 ('https://pypi.python.org/pypi/scikit-learn/0.16.1'))
			return

	except ImportError:
		sys.stderr.write('scikit-learn package not installed: install version 0.16.1 from here: %s \n' %
						 ('https://pypi.python.org/pypi/scikit-learn/0.16.1'))
		return

	#Rule Set 2 scoring models
	preinstalled_V3_full_model,preinstalled_V3_nopos_model = 'Rule_Set_2_scoring/saved_models/V3_model_full.pickle',\
															 'Rule_Set_2_scoring/saved_models/V3_model_nopos.pickle'
	V3_model_full = resource_filename(__name__, preinstalled_V3_full_model)
	V3_model_nopos = resource_filename(__name__, preinstalled_V3_nopos_model)

	if indir2 and filename2:
		db1 = indir + '/' + filename
		db2 = indir2 + '/' + filename2
		database_path_verification(db1)
		database_path_verification(db2)
		database_check(db1)
		database_check(db2)
		sys.stdout.write('database to database cutting efficiency score transfer functions called \n')
		ces_database_to_database_transfer.core(indir,filename,indir2,filename2)
		os.rename('newdatabase_w_transfered_scores.bam', 'newdatabase_w_transfered_efficiency_scores.bam')
		os.rename('newdatabase_w_transfered_scores.bam.bai', 'newdatabase_w_transfered_efficiency_scores.bam.bai')
		sys.stdout.write('finished \n')
	else:
		database_path_verification(indir + '/' +filename)
		database_check(indir + '/' +filename)
		sys.stdout.write('Rule Set 2 cutting efficiency score computation and insertion functions called \n')
		ces_compute.core(indir, filename, fasta, aa_cut, per_peptide, V3_model_nopos, V3_model_full)
		score = 'CutEffScore_' + filename.replace('.bam','.txt')
		ces_insert.core(indir,filename,score)
		sys.stdout.write('finished \n')

if __name__ == '__main__':
	main()


