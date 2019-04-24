__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""computes and inserts CFD and/or Hsu cutting specificity scores into Cas9 sgRNA database"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import os
import argparse

from pkg_resources import resource_exists, resource_filename

import ces_insert
import ces_database_to_database_transfer
import tss_compute

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def kmers_counted_file_trawl(search_directory):
    """get or generate kmer counts file for specificity scoring

    :param search_directory: absolute filepath to directory with file(s) having _all_kmers.txt.gz or _all_kmers_counted.txt suffix
    :return: filepath to _all_kmers_counted.txt file

    """
    kmers_lst, kmers_counted_lst = kmer_file_trawler(search_directory)
    if kmers_counted_lst:
        if len(kmers_counted_lst) == 1:
            kmer_counts_file = kmers_counted_lst[0]
            sys.stdout.write('_kmers_counted.txt file found\n')
            return kmer_counts_file
        else:
            sys.stderr.write('ERROR: %s many files with %s suffix, only one can be present in directory\n' % (
            len(kmers_counted_lst), '_kmers_counted.txt'))
            return
    elif kmers_lst:
        if len(kmers_lst) == 1:
            sys.stdout.write('generating %s_counted.txt\n' % (kmers_lst[0].split('.')[0]))
            cmd = "gzip -cd %s | awk '{print $1}' - | sort - | uniq -c - > %s_counted.txt" % (
                kmers_lst[0], kmers_lst[0].split('.')[0])
            os.system(cmd)
            kmer_counts_file = ('%s_counted.txt' % (kmers_lst[0].split('.')[0]))
            return kmer_counts_file
        else:
            sys.stderr.write('ERROR: %s many files with %s suffix, only one can be present in directory\n' % (
                len(kmers_lst), '_all_kmers.txt.gz'))
            return
    else:
        sys.stderr.write(
            'ERROR: Neither X_all_kmers.txt.gz nor X_all_kmers_counted.txt file not located in %s, specificity scoring not engaged\n' % (
            search_directory))
        return

def kmer_file_trawler(file_dir):
    """finds all .txt files in a given directory

    Input:
    file_dir: filepath to query directory

    """
    kmers_lst,kmers_counted_lst = [],[]
    for file in os.listdir(file_dir):
        if file.endswith('_kmers.txt.gz'):
            kmers_lst.append('%s/%s' % (file_dir,file))
        elif file.endswith('_kmers_counted.txt'):
            kmers_counted_lst.append('%s/%s' % (file_dir,file))

    return kmers_lst,kmers_counted_lst

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--database_directory', help='absolute filepath to directory hosting sgRNA database with '
														   'no cutting specificity scores', required=True)
	parser.add_argument('-n', '--database_name', help='filename of sgRNA database with no cutting specificity scores',
						required=True)
	parser.add_argument('-k', '--kmers_file',help='absolute filepath to directory hosting either X_all_kmers.txt.gz or X_all_kmers_counted.txt')
	parser.add_argument('-f', '--fasta_filepath',
						help='absolute filepath to organism FASTA file. Index for FASTA should '
							 'also be present in same directory (.fai). Must be single aggregate FASTA file for organism')
	parser.add_argument('-d2', '--database_directory2',
						help='absolute filepath to directory hosting sgRNA database with cutting '
							 'specificity scores already included and generated with all the '
							 'same run parameters as another database except for -d. This '
							 'parameter, coupled with -n2, allows for the transfer of '
							 'cutting efficiency scores between two sgRNA databases that '
							 'differ only in how many mismatches off-targets are '
							 'enumerated to')
	parser.add_argument('-n2', '--database_name2',
						help='filename of sgRNA database with cutting specificity scores already'
							 ' included and generated with all the same run parameters as '
							 'another database except for -d. This parameter, coupled with -d2,'
							 ' allows for the transfer of cutting efficiency scores between two'
							 ' sgRNA databases that differ only in how many mismatches '
							 'off-targets are enumerated to')

	args = parser.parse_args()
	indir = args.database_directory
	filename = args.database_name
	kmers = args.kmers_file
	fasta = args.fasta_filepath
	indir2 = args.database_directory2
	filename2 = args.database_name2

	return indir, filename, kmers, fasta, indir2, filename2

def database_check(database):
	"""function that verifies that a sgRNA database meets the requirements for CFD and/or Hsu cutting specificity
		score computation

		Input:
		database: output of GuideScan.processor. This will be the BAM file database of sgRNAs that meet the user's
		predefined parameters. The check is done by ensuring proper run parameters for CFD scoring are present in
		bam header.

		Note:
		For a database to be eligible for CFD cutting specificity scoring, it must be composed of sgRNAs
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
			sys.stderr.write('ERROR: database %s does not meet requirements for CFD cutting specificity scoring \n'
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
			sys.stdout.write('MSG: database %s meets requirements for CFD cutting specificity score computation \n'
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

	#user inputs
	indir, filename, kmers, fasta, indir2, filename2 = arg_parser()

	#CFD scoring pickles
	mm_score_pickle, pam_score_pickle = 'CFD_scoring/mismatch_score.pkl','CFD_scoring/pam_scores.pkl'
	mismatch_score = resource_filename(__name__, mm_score_pickle)
	pam_score = resource_filename(__name__, pam_score_pickle)

	#cutting specificity scoring
	if indir2 and filename2:
		db1 = indir + '/' + filename
		db2 = indir2 + '/' + filename2
		database_path_verification(db1)
		database_path_verification(db2)
		database_check(db1)
		database_check(db2)
		sys.stdout.write('database to database cutting efficiency score transfer functions called \n')
		ces_database_to_database_transfer.core(indir, filename, indir2, filename2)
		os.rename('newdatabase_w_transfered_scores.bam','newdatabase_w_transfered_specificity_scores.bam')
		os.rename('newdatabase_w_transfered_scores.bam.bai','newdatabase_w_transfered_specificity_scores.bam.bai')
		sys.stdout.write('finished \n')
	else:
		# kmers file
		kmer_counts_file = kmers_counted_file_trawl(kmers)
		if kmer_counts_file:
			database_path_verification(indir + '/' + filename)
			database_check(indir + '/' + filename)
			sys.stdout.write('CFD specificity score computation and insertion functions called \n')
			tss_compute.core(indir,filename,fasta,mismatch_score,pam_score,kmer_counts_file)
			score = 'scoring_outfile_' + filename.replace('.bam', '.txt')
			ces_insert.core(indir, filename, score)
			os.rename('cutting_efficiency_scores_added_' + filename, 'cutting_specificity_scores_added_' + filename)
			os.rename('cutting_efficiency_scores_added_' + filename + '.bai', 'cutting_specificity_scores_added_' + filename + '.bai')
			sys.stdout.write('finished \n')
		else:
			sys.stderr.write('ERROR: specificity scores not computed\n')
			return

if __name__ == '__main__':
	main()