__author__ = 'Alexendar Perez & Yuri Pritykin & Sagar Chhanagawala'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""pickles gRNA sets from GuideScan gRNA databases"""

#################
#               #
#   Libraries   #
#               #
#################

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

class gRNA_Chained_Hash_Table:
	"""class which reads in GuideScan gRNA file where first field is gRNA target sequence and returns chained hash table"""

	def __init__(self,infile,capacity):
		self.infile = infile
		self.hash_table = None
		self.capacity = capacity
		self.slots = []
		for i in range(capacity):
			self.slots.append([])

	def __len__(self):
		count = 0
		for i in self.slots:
			count += len(i)
		return count

	def __str__(self):
		content = ''
		for i in self.slots:
			'%s' % str(i)
		return content

	def __getitem__(self, item):
		slot = sum((map(hash, item))) % self.capacity
		if item in self.slots[slot]:
			return 0
		else:
			return -1

	def _insert(self,key):
		#print key,map(hash,key)
		slot = sum((map(hash,key))) % self.capacity
		if key in self.slots[slot]:
			return -1
		else:
			self.slots[slot].append(key)
			return slot

	def _make_hash_table(self):
		hash_table = None
		with open(self.infile,'rb') as file_in:
			for line in file_in:
				clean_line = line.lstrip().rstrip()
				parts = clean_line.split()
				hash_table = self._insert(parts[0])

		self.hash_table = hash_table
		sys.stdout.write('chained hash table constructed\n')

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--infile',help='absolute filepath to input file',required=True)
	parser.add_argument('-o','--outdir',help='absolute filepath to output directory',required=True)
	parser.add_argument('-id','--organism_id',help='organism assembly id, (ie: human = hg38)',required=True)

	args = parser.parse_args()
	gRNA_file = args.infile
	outdir = args.outdir
	organism = args.organism_id

	return gRNA_file,outdir,organism

def construct_and_serialize_gRNA_sets(gRNA_file,outdir,organism='model_organism'):
	"""construct and serialize a gRNA database for inclusion queries

	Inputs:
	gRNA_file: absolute filepath to user inputted file containing gRNAs to be queried for inclusion in GuideScan database
	outdir: absolute filepath to output directory for writeout

	Outputs:
	a python pickle file which contains the set of gRNA target sequences present for in a GuideScan database for a given
	organism

	"""
	#create set
	beta = gRNA_set(gRNA_file)
	beta._make_set()
	#serialize set
	with open('%s/gRNA_set_%s.pkl' % (outdir,organism), 'wb') as outfile:
		cPickle.dump(beta, outfile)

	sys.stdout.write('%s written for %s\n' % ('%s/gRNA_set_%s.pkl' % (outdir,organism),gRNA_file))

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

	# user inputs
	gRNA_file, outdir, organism = arg_parser()

	# construct pickle
	construct_and_serialize_gRNA_sets(gRNA_file, outdir, organism)

	sys.stdout.write('pickle file for %s generated\n' % (organism))

if __name__ == '__main__':
	main()