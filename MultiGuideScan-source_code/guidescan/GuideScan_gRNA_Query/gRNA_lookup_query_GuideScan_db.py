__author__ = 'Alexendar Perez & Yuri Pritykin & Sagar Chhanagawala'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""allow for user to query their gRNA sequence for inclusion inside GuideScan database"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import re
import sys
import cPickle
import argparse

#############
#           #
#   Class   #
#           #
#############

class gRNA_set:
	"""class which reads in GuideScan gRNA file where first field is gRNA target sequence and returns set object"""

	def __init__(self,infile):
		self.f = infile
		self.set = set()

	def __getitem__(self, item):
		if item in self.set:
			return 0
		else:
			return 1

	def _insert(self,item):
		return self.set.add(item)

	def _make_set(self):
		with open(self.f,'rb') as file_in:
			for line in file_in:
				clean_line = line.lstrip().rstrip()
				parts = clean_line.split()
				self._insert(parts[0])

		sys.stdout.write('gRNA set contructed\n')

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pickles_dir',help='absolute filepath to directory hosting gRNA pickle files',required=True)
	parser.add_argument('-i','--infile',help='absolute filepath to input file with queried gRNA sequences',required=True)
	parser.add_argument('--organism',help='string object detailing which organism under query: dm6,mm10,hg38,danRer10,ce11,SacCerv',required=True)
	parser.add_argument('-o','--outdir',help='absolute filepath to output directory for writeout',required=True)
	parser.add_argument('--selective',help='if not None then only deserialize queried organism pickle; default is None',default='')

	args = parser.parse_args()
	pickles = args.pickles_dir
	infile = args.infile
	organism = args.organism
	outdir = args.outdir
	selective = args.selective

	return pickles,infile,organism,outdir,selective

def deserialize_pickles(outdir,organism,selective):
	"""deserialize gRNA pickles

	:param outdir: absolute filepath to directory hosting pickles
	:param organism: string identifier of organism (ie: hg38)
	:param selective: if not None then returns first pickle object
	:return: dictionary object with name of organism as key and deserialized pickle sets of gRNAs as values

	Convention:
	selective is for individual queries, when set to None is ideal for implementation into web-interface

	"""

	# collect pickle files
	file_lst = file_trawler(outdir, '.pkl', 'gRNA_set')

	# deserialize one or all pickles
	if selective:
		query_file = '%s/%s_%s%s' % (outdir, 'gRNA_set', organism, '.pkl')
		if query_file in file_lst:
			file_lst = [query_file]
		else:
			sys.stderr.write('%s did not correspond to a saved pickled file, deserialize all pickles\n' % (query_file ))
			selective = None

	# models list
	model_dict = {}

	# deserialize pickles
	for item in file_lst:
		org_name = item.split('/')[-1].split('_')[-1].rstrip('.pkl')
		if org_name == 'dm6':
			with open(item, 'rb') as infile:
				dm6 = cPickle.load(infile)
				if selective:
					model_dict['dm6'] = dm6
					sys.stdout.write('dm6 pickle deserialzed\n')
					return model_dict
				else:
					model_dict['dm6'] = dm6
					sys.stdout.write('dm6 pickle deserialzed\n')
		elif org_name == 'mm10':
			with open(item, 'rb') as infile:
				mm10 = cPickle.load(infile)
				if selective:
					model_dict['mm10'] = mm10
					sys.stdout.write('mm10 pickle deserialzed\n')
					return model_dict
				else:
					model_dict['mm10'] = mm10
					sys.stdout.write('mm10 pickle deserialzed\n')
		elif org_name == 'hg38':
			with open(item, 'rb') as infile:
				hg38 = cPickle.load(infile)
				if selective:
					model_dict['hg38'] = hg38
					sys.stdout.write('hg38 pickle deserialzed\n')
					return model_dict
				else:
					model_dict['hg38'] = hg38
					sys.stdout.write('hg38 pickle deserialzed\n')
		elif org_name == 'ce11':
			with open(item, 'rb') as infile:
				ce11 = cPickle.load(infile)
				if selective:
					model_dict['ce11'] = ce11
					sys.stdout.write('ce11 pickle deserialzed\n')
					return model_dict
				else:
					model_dict['ce11'] = ce11
					sys.stdout.write('ce11 pickle deserialzed\n')
		elif org_name == 'danRer10':
			with open(item, 'rb') as infile:
				danRer10 = cPickle.load(infile)
				if selective:
					model_dict['danRer10'] = danRer10
					sys.stdout.write('danRer10 pickle deserialzed\n')
					return model_dict
				else:
					model_dict['danRer10'] = danRer10
					sys.stdout.write('danRer10 pickle deserialzed\n')
		elif org_name == 'SacCerv':
			with open(item, 'rb') as infile:
				SacCerv = cPickle.load(infile)
				if selective:
					model_dict['SacCerv'] = SacCerv
					sys.stdout.write('SacCerv pickle deserialzed\n')
					return model_dict
				else:
					model_dict['SacCerv'] = SacCerv
					sys.stdout.write('SacCerv pickle deserialzed\n')
		else:
			sys.stdout.write('unknown organism pickle: %s\n' % (org_name))
			with open(item, 'rb') as infile:
				model = cPickle.load(infile)
				if selective:
					model_dict['X'] = model
					return model_dict
				else:
					model_dict['X'] = model

	sys.stdout.write('deserialization complete\n')
	return model_dict

def file_trawler(indir, suffix, conserved_string=''):
	"""locates all files with suffix ending in specified directory

	Inputs:
	indir = absolute filepath to directory where files with desired suffix are located
	suffix = file suffix of desired files (ie: .txt)
	conserved_string: conserved string corresponding to files to consider; default is none

	Output:
	file_lst: list object with the filepaths that conform to the suffix and conserved_string requirements

	"""
	files = os.listdir(indir)
	file_lst = []
	for file in files:
		if file.endswith(suffix):
			if conserved_string:
				if conserved_string in file:
					file_lst.append('%s/%s' % (indir, file))
				else:
					continue
			else:
				file_lst.append('%s/%s' % (indir, file))
		else:
			continue

	sys.stdout.write('files ending with %s and conserved string "%s" were appended to file list \n' % (suffix,conserved_string))
	return file_lst

def database_inclusion_query_and_writeout(gamma,gRNA_in_file,outdir,query_organism):
	"""determines if gRNAs in user inputted file are contained within GuideScan database

		Inputs:
		gamma: output of gRNA_set class two set command: 1.) gamma = gRNA_set(gRNA_from_GuideScan_file) 2.) gamma._make_set()
		gRNA_in_file: absolute filepath to user inputted file containing gRNAs to be queried for inclusion in GuideScan database
		outdir: absolute filepath to output directory for writeout

		Outputs:
		A text file is created in outdir which is composed of two fields. The first field contains the the queired gRNA
		while the second field includes a string indicating where the gRNA was found in the GuideScan database or was not

		"""
	with open('%s/GuideScan_gRNA_inclusion_%s.txt' % (outdir,query_organism), 'wb') as outfile:
		with open(gRNA_in_file, 'rb') as infile:
			for line in infile:
				if line != '\n':
					clean_line = line.lstrip().rstrip()
					parts = clean_line.split()
					if len(parts[0]) == 20:
						if re.match("^[A-Z]*$",str(parts[0]).upper()):
							if gamma.__getitem__(parts[0]) == 0:
								outfile.write('%s\tIncluded in GuideScan database\n' % (parts[0]))
							else:
								outfile.write('%s\tExcluded from GuideScan database\n' % (parts[0]))
						else:
							sys.stdout.write('pure character string not encountered with %s, skipping\n' % (parts[0]))
							continue
					else:
						sys.stdout.write('gRNA complementary sequence is not 20bp, database is of 20bp gRNAs\n')
						continue
				else:
					sys.stdout.write('newline character encountered, skipping\n')
					continue

	sys.stdout.write('GuideScan inclusion file written to %s\n' % '%s/GuideScan_gRNA_inclusion.txt' % (outdir))

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

	"""
	pickles = '/Users/pereza1/Projects/Ventura/CRISPR/data/sorted_gRNA_text_files/pickles'
	infile = '/Users/pereza1/Projects/Ventura/CRISPR/data/sorted_gRNA_text_files/dm6_test/positive_samples.txt'
	organism = 'dm6'
	outdir = '/Users/pereza1/Desktop'
	"""

	pickles,infile,organism,outdir,selective = arg_parser()


	organism_dictionary = {'hg38': 'hg38', 'human': 'hg38', 'Hg38': 'hg38', 'HG38': 'hg38',
						   'mm10': 'mm10', 'mouse': 'mm10', 'Mm10': 'mm10', 'MM10': 'mm10',
						   'dm6': 'dm6', 'fly': 'dm6', 'drosophila': 'dm6', 'fruit_fly': 'dm6', 'Dm6': 'dm6', 'DM6': 'dm6',
						   'ce11': 'ce11', 'worm': 'ce11', 'c.eleigan': 'ce11', 'Ce11': 'ce11', 'cell': 'ce11',
						   'CE11': 'cd11',
						   'danRer10': 'danRer10', 'zebrafish': 'danRer10', 'danrer10': 'danRer10', 'DANRER10': 'danRer10',
						   'SacCerv': 'SacCerv', 'yeast': 'SacCerv', 'saccerv': 'SacCerv', 'SACCERV': 'SacCerv'}

	if organism_dictionary.has_key(organism):
		sys.stdout.write('%s called %s GuideScan database\n' % (organism,organism_dictionary[organism]))
		organism = organism_dictionary[organism]
	else:
		sys.stderr.write('WARNING: %s is not a recognized species name\n' % (organism))
		return

	gRNAs = deserialize_pickles(pickles,organism,selective)
	if gRNAs:
		database_inclusion_query_and_writeout(gRNAs[organism],infile,outdir,organism)
		return
	else:
		sys.stderr.write('WARNING: pickles not deseralized, ensure proper filepath to pickles\n')
		return

if __name__ == '__main__':
	main()



