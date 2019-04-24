__author__ = 'Alex Perez and Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Computes cutting efficiency scores of sgRNAs

Note: Cutting efficiency is computed using the Rule 2 method devised by Doench et al. Nature Biotechnology 2016.
Functions under --Rule 2 Functions-- are used wholesale and without modification from original paper. Functions under
--Auxillary Functions-- are written by authors. Author written functions are summarized below:

#####################################
#                                   #
#   Auxillary Functions summaries   #
#                                   #
#####################################

serial_on_target_doench(seq,aa_cut,per_peptide,model) = DSF, ensures GG is located in the sterotyped NGG positions for
														the PAM sequences of the 30mer
donech_scoring_parameters(aa_cut,per_peptide,V3_model_nopos,V3_model_full) = Load the correct Donech scoring model
bamfile_info_extraction(indir,filename) = Extract chromosome, coordinate, strand, original sequence from BAM file
fasta(fasta_file) = Load organism fasta file for use in pyfaidx module
cutting_efficiency_w_doench_rule2(outfile,indir,filename,org,aa_cut,per_peptide,model) = Compute cutting efficiency
																						 using Rule 2 model
"""

#################
#               #
#   Libraries   #
#               #
#################

from pyfaidx import Fasta
from Bio import Seq
import os
import argparse
import pickle
import pandas
import time
import numpy as np
import Bio.SeqUtils as SeqUtil
import sys
import Bio.SeqUtils.MeltingTemp as Tm
import itertools
import sklearn

#############################
#                           #
#      Rule 2 Functions     #
#                           #
#############################

#these functions are from Doench et al. Nature Biotechnology 2016

def get_gene_sequence(gene_name):
	try:
		gene_file = '../gene_sequences/%s_sequence.txt' % gene_name
		with open(gene_file, 'rb') as f:
			seq = f.read()
			seq = seq.replace('\r\n', '')
	except:
		raise Exception("could not find gene sequence file %s, please see examples and generate one for your gene as needed, with this filename" % gene_file)

	return seq

def concatenate_feature_sets(feature_sets):
	'''
	Given a dictionary of sets of features, each in a Pandas.DataFrame,
	concatenate them together to form one big np.array, and get the dimension
	of each set
	Returns: inputs, dim
	'''
	assert feature_sets != {}, "no feature sets present"
	F = feature_sets[feature_sets.keys()[0]].shape[0]
	for set in feature_sets.keys():
		F2 = feature_sets[set].shape[0]
		assert F == F2, "not same # individuals for features %s and %s" % (feature_sets.keys()[0], set)

	N = feature_sets[feature_sets.keys()[0]].shape[0]
	inputs = np.zeros((N, 0))
	feature_names = []
	dim = {}
	dimsum = 0
	for set in feature_sets.keys():
		inputs_set = feature_sets[set].values
		dim[set] = inputs_set.shape[1]
		dimsum += dim[set]
		inputs = np.hstack((inputs, inputs_set))
		feature_names.extend(feature_sets[set].columns.tolist())

	return inputs, dim, dimsum, feature_names

def predict(seq, aa_cut=0, percent_peptide=0, model=None, model_file=None):
	assert not (model is None and model_file is None), "you have to specify either a model or a model_file"

	if model is None:
		try:
			with open(model_file, 'rb') as f:
				model, learn_options = pickle.load(f)
		except:
			raise Exception("could not find model stored to file %s" % model_file)
	else:
		model, learn_options = model

	learn_options["V"] = 2

	# Y, feature_sets, target_genes, learn_options, num_proc = setup(test=False, order=2, learn_options=learn_options, data_file=test_filename)
	# inputs, dim, dimsum, feature_names = pd.concatenate_feature_sets(feature_sets)

	Xdf = pandas.DataFrame(columns=[u'30mer', u'Strand'], data=[[seq, 'NA']])
	gene_position = pandas.DataFrame(columns=[u'Percent Peptide', u'Amino Acid Cut position'], data=[[percent_peptide, aa_cut]])
	feature_sets = featurize_data(Xdf, learn_options, pandas.DataFrame(), gene_position)
	inputs, dim, dimsum, feature_names = concatenate_feature_sets(feature_sets)

	# call to scikit-learn, returns a vector of predicted values
	return model.predict(inputs)[0]

def featurize_data(data, learn_options, Y, gene_position):
	'''
	assumes that data contains the 30mer
	returns set of features from which one can make a kernel for each one
	'''

	#print "Constructing features..."
	t0 = time.time()

	feature_sets = {}

	if learn_options["nuc_features"]:
		# spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
		get_all_order_nuc_features(data['30mer'], feature_sets, learn_options, learn_options["order"], max_index_to_use=30)

	if learn_options["gc_features"]:
		gc_above_10, gc_below_10, gc_count = gc_features(data)
		feature_sets['gc_above_10'] = pandas.DataFrame(gc_above_10)
		feature_sets['gc_below_10'] = pandas.DataFrame(gc_below_10)
		feature_sets['gc_count'] = pandas.DataFrame(gc_count)

	if learn_options["include_gene_position"]:
		# gene_position_columns = ["Amino Acid Cut position", "Percent Peptide", "Nucleotide cut position"]
		# gene_position_columns = ["Percent Peptide", "Nucleotide cut position"]

		for set in gene_position.columns:
			set_name = set
			feature_sets[set_name] = pandas.DataFrame(gene_position[set])
		feature_sets["Percent Peptide <50%"] = feature_sets["Percent Peptide"] < 50

	if learn_options["include_gene_effect"]:
		print "including gene effect"
		gene_names = Y['Target gene']
		enc = sklearn.preprocessing.OneHotEncoder()
		label_encoder = sklearn.preprocessing.LabelEncoder()
		label_encoder.fit(gene_names)
		one_hot_genes = np.array(enc.fit_transform(label_encoder.transform(gene_names)[:, None]).todense())
		feature_sets["gene effect"] = pandas.DataFrame(one_hot_genes,
													   columns=["gene_%d" % i for i in range(one_hot_genes.shape[1])], index=gene_names.index)

	if learn_options['include_known_pairs']:
		feature_sets['known pairs'] = pandas.DataFrame(Y['test'])

	if learn_options["include_NGGX_interaction"]:
		feature_sets["NGGX"] = NGGX_interaction_feature(data)

	if learn_options["include_Tm"]:
		feature_sets["Tm"] = Tm_feature(data)

	if learn_options["include_sgRNAscore"]:
		feature_sets["sgRNA Score"] = pandas.DataFrame(data["sgRNA Score"])

	if learn_options["include_drug"]:
		# feature_sets["drug"] = pandas.DataFrame(data["drug"])
		drug_names = Y.index.get_level_values('drug').tolist()
		enc = sklearn.preprocessing.OneHotEncoder()
		label_encoder = sklearn.preprocessing.LabelEncoder()
		label_encoder.fit(drug_names)
		one_hot_drugs = np.array(enc.fit_transform(label_encoder.transform(drug_names)[:, None]).todense())
		feature_sets["drug"] = pandas.DataFrame(one_hot_drugs, columns=["drug_%d" % i for i in range(one_hot_drugs.shape[1])], index=drug_names)

	if learn_options['include_strand']:
		feature_sets['Strand effect'] = (pandas.DataFrame(data['Strand']) == 'sense')*1

	if learn_options["include_gene_feature"]:
		feature_sets["gene features"] = gene_feature(Y, data, learn_options)

	if learn_options["include_gene_guide_feature"] > 0:
		tmp_feature_sets = gene_guide_feature(Y, data, learn_options)
		for key in tmp_feature_sets:
			feature_sets[key] = tmp_feature_sets[key]

	if learn_options["include_microhomology"]:
		feature_sets["microhomology"] = get_micro_homology_features(Y['Target gene'], learn_options, data)

	t1 = time.time()
	#print "\t\tElapsed time for constructing features is %.2f seconds" % (t1-t0)

	check_feature_set_dimensions(feature_sets)

	if learn_options['normalize_features']:
		feature_sets = normalize_feature_sets(feature_sets)
		check_feature_set_dimensions(feature_sets)

	return feature_sets

def check_feature_set_dimensions(feature_sets):
	'''
	Ensure the # of people is the same in each feature set
	'''
	N = None
	for ft in feature_sets.keys():
		N2 = feature_sets[ft].shape[0]
		if N is None:
			N = N2
		else:
			assert N == N2, "# of individuals do not match up across feature sets"

def NGGX_interaction_feature(data):
	'''
	assuming 30-mer, grab the NGGX _ _ positions, and make a one-hot
	encoding of the NX nucleotides yielding 4x4=16 features
	'''
	sequence = data['30mer'].values
	feat_NX = pandas.DataFrame()
	# check that GG is where we think
	for seq in sequence:
		if seq[25:27] != "GG":
			raise Exception("expected GG but found %s" % seq[25:27])
		NX = seq[24]+seq[27]
		NX_onehot = nucleotide_features(NX,order=2, include_pos_independent=False, max_index_to_use=2, prefix="NGGX")
		# NX_onehot[:] = np.random.rand(NX_onehot.shape[0]) ##TESTING RANDOM FEATURE
		feat_NX = pandas.concat([feat_NX, NX_onehot], axis=1)
	return feat_NX.T

def get_all_order_nuc_features(data, feature_sets, learn_options, maxorder, max_index_to_use, prefix=""):
	for order in range(1, maxorder+1):
		#print "\t\tconstructing order %s features" % order
		nuc_features_pd, nuc_features_pi = apply_nucleotide_features(data, order, learn_options["num_proc"],
																	 include_pos_independent=True, max_index_to_use=max_index_to_use, prefix=prefix)
		feature_sets['%s_nuc_pd_Order%i' % (prefix, order)] = nuc_features_pd
		if learn_options['include_pi_nuc_feat']:
			feature_sets['%s_nuc_pi_Order%i' % (prefix, order)] = nuc_features_pi
		#print "\t\t\t\t\t\t\tdone"

def countGC(s):
	'''
	GC content for only the 20mer, as per the Doench paper/code
	'''
	assert len(s) == 30, "seems to assume 30mer"
	return len(s[5:25].replace('A', '').replace('T', ''))

def SeqUtilFeatures(data):
	'''
	assuming '30-mer'is a key
	get melting temperature features from:
		0-the 30-mer ("global Tm")
		1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
		2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
		3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
	'''
	sequence = data['30mer'].values
	num_features = 1
	featarray = np.ones((sequence.shape[0], num_features))
	for i, seq in enumerate(sequence):
		assert len(seq) == 30, "seems to assume 30mer"
		featarray[i, 0] = SeqUtil.molecular_weight(str(seq))

	feat = pandas.DataFrame(pandas.DataFrame(featarray))
	return feat

def organism_feature(data):
	'''
	Human vs. mouse
	'''
	organism = np.array(data['Organism'].values)
	feat = pandas.DataFrame(pandas.DataFrame(featarray))
	import ipdb; ipdb.set_trace()
	return feat

def get_micro_homology_features(gene_names, learn_options, X):
	# originally was flipping the guide itself as necessary, but now flipping the gene instead

	print "building microhomology features"
	feat = pandas.DataFrame(index=X.index)
	feat["mh_score"] = ""
	feat["oof_score"] = ""

	#with open(r"tmp\V%s_gene_mismatches.csv" % learn_options["V"],'wb') as f:
	if True:
		# number of nulceotides to take to the left and right of the guide
		k_mer_length_left = 9
		k_mer_length_right = 21
		for gene in gene_names.unique():
			gene_seq = Seq.Seq(get_gene_sequence(gene)).reverse_complement()
			guide_inds = np.where(gene_names.values == gene)[0]
			print "getting microhomology for all %d guides in gene %s" % (len(guide_inds), gene)
			for j, ps in enumerate(guide_inds):
				guide_seq = Seq.Seq(X['30mer'][ps])
				strand = X['Strand'][ps]
				if strand=='sense':
					gene_seq = gene_seq.reverse_complement()
				# figure out the sequence to the left and right of this guide, in the gene
				ind = gene_seq.find(guide_seq)
				if ind==-1:
					gene_seq = gene_seq.reverse_complement()
					ind = gene_seq.find(guide_seq)
					#assert ind != -1, "still didn't work"
					#print "shouldn't get here"
				else:
					#print "all good"
					pass
				#assert ind != -1, "could not find guide in gene"
				if ind==-1:
					#print "***could not find guide %s for gene %s" % (str(guide_seq), str(gene))
					#if.write(str(gene) + "," + str(guide_seq))
					mh_score = 0
					oof_score = 0
				else:
					#print "worked"

					assert gene_seq[ind:(ind+len(guide_seq))]==guide_seq, "match not right"
					left_win = gene_seq[(ind - k_mer_length_left):ind]
					right_win = gene_seq[(ind + len(guide_seq)):(ind + len(guide_seq) + k_mer_length_right)]

					#if strand=='antisense':
					#    # it's arbitrary which of sense and anti-sense we flip, we just want
					#    # to keep them in the same relative alphabet/direction
					#    left_win = left_win.reverse_complement()
					#    right_win = right_win.reverse_complement()
					assert len(left_win.tostring())==k_mer_length_left
					assert len(right_win.tostring())==k_mer_length_right

					sixtymer = str(left_win) + str(guide_seq) + str(right_win)
					assert len(sixtymer)==60, "should be of length 60"
					mh_score, oof_score = microhomology.compute_score(sixtymer)

				feat.ix[ps,"mh_score"] = mh_score
				feat.ix[ps,"oof_score"] = oof_score
			print "computed microhomology of %s" % (str(gene))

	return pandas.DataFrame(feat, dtype='float')

def local_gene_seq_features(gene_names, learn_options, X):

	print "building local gene sequence features"
	feat = pandas.DataFrame(index=X.index)
	feat["gene_left_win"] = ""
	feat["gene_right_win"] = ""

	# number of nulceotides to take to the left and right of the guide
	k_mer_length = learn_options['include_gene_guide_feature']
	for gene in gene_names.unique():
		gene_seq = Seq.Seq(get_gene_sequence(gene)).reverse_complement()
		for ps in np.where(gene_names.values==gene)[0]:
			guide_seq = Seq.Seq(X['30mer'][ps])
			strand = X['Strand'][ps]
			if strand=='sense':
				guide_seq = guide_seq.reverse_complement()
				#gene_seq = gene_seq.reverse_complement()
			# figure out the sequence to the left and right of this guide, in the gene
			ind = gene_seq.find(guide_seq)
			if ind ==-1:
				#gene_seq = gene_seq.reverse_complement()
				#ind = gene_seq.find(guide_seq)
				assert ind != -1, "could not find guide in gene"
			assert gene_seq[ind:(ind+len(guide_seq))]==guide_seq, "match not right"
			left_win = gene_seq[(ind - k_mer_length):ind]
			right_win = gene_seq[(ind + len(guide_seq)):(ind + len(guide_seq) + k_mer_length)]

			if strand=='antisense':
				# it's arbitrary which of sense and anti-sense we flip, we just want
				# to keep them in the same relative alphabet/direction
				left_win = left_win.reverse_complement()
				right_win = right_win.reverse_complement()
			assert not left_win.tostring()=="", "k_mer_context, %s, is too large" % k_mer_length
			assert not left_win.tostring()=="", "k_mer_context, %s, is too large" % k_mer_length
			assert len(left_win)==len(right_win), "k_mer_context, %s, is too large" % k_mer_length
			feat.ix[ps,"gene_left_win"] = left_win.tostring()
			feat.ix[ps,"gene_right_win"] = right_win.tostring()
		print "featurizing local context of %s" % (gene)

	feature_sets = {}
	get_all_order_nuc_features(feat["gene_left_win"], feature_sets, learn_options, learn_options["order"], max_index_to_use=sys.maxint, prefix="gene_left_win")
	get_all_order_nuc_features(feat["gene_right_win"], feature_sets, learn_options, learn_options["order"], max_index_to_use=sys.maxint, prefix="gene_right_win")
	return feature_sets

def gene_feature(Y, X, learn_options):
	'''
	Things like the sequence of the gene, the DNA Tm of the gene, etc.
	'''

	gene_names = Y['Target gene']

	gene_length = np.zeros((gene_names.values.shape[0], 1))
	gc_content = np.zeros((gene_names.shape[0], 1))
	temperature = np.zeros((gene_names.shape[0], 1))
	molecular_weight = np.zeros((gene_names.shape[0], 1))

	for gene in gene_names.unique():
		seq = get_gene_sequence(gene)
		gene_length[gene_names.values==gene] = len(seq)
		gc_content[gene_names.values==gene] = SeqUtil.GC(seq)
		temperature[gene_names.values==gene] = Tm.Tm_staluc(seq, rna=False)
		molecular_weight[gene_names.values==gene] = SeqUtil.molecular_weight(seq, 'DNA')

	all = np.concatenate((gene_length, gc_content, temperature, molecular_weight), axis=1)
	df = pandas.DataFrame(data=all, index=gene_names.index, columns=['gene length',
																	 'gene GC content',
																	 'gene temperature',
																	 'gene molecular weight'])
	return df

def gene_guide_feature(Y, X, learn_options):
	#features, which are related to parts of the gene-local to the guide, and
	#possibly incorporating the guide or interactions with it

	#expensive, so pickle if necessary
	gene_file = r"..\data\gene_seq_feat_V%s_km%s.ord%s.pickle" % (learn_options['V'], learn_options['include_gene_guide_feature'], learn_options['order'])

	if False: #os.path.isfile(gene_file): #while debugging, comment out
		print "loading local gene seq feats from file %s" % gene_file
		with open(gene_file, "rb") as f: feature_sets = pickle.load(f)
	else:
		feature_sets = local_gene_seq_features(Y['Target gene'], learn_options, X)
		print "writing local gene seq feats to file %s" % gene_file
		with open(gene_file, "wb") as f: pickle.dump(feature_sets, f)

	return feature_sets

def gc_cont(seq):
	return (seq.count('G') + seq.count('C'))/float(len(seq))

def Tm_feature(data):
	'''
	assuming '30-mer'is a key
	get melting temperature features from:
		0-the 30-mer ("global Tm")
		1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
		2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
		3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
	'''
	sequence = data['30mer'].values
	featarray = np.ones((sequence.shape[0],4))
	for i, seq in enumerate(sequence):
		if seq[25:27]!="GG":
			raise Exception("expected GG but found %s" % seq[25:27])
		rna = False
		featarray[i,0] = Tm.Tm_staluc(seq, rna=rna)        #30mer Tm
		featarray[i,1] = Tm.Tm_staluc(seq[20:25], rna=rna) #5nts immediately proximal of the NGG PAM
		featarray[i,2] = Tm.Tm_staluc(seq[12:20], rna=rna)   #8-mer
		featarray[i,3] = Tm.Tm_staluc(seq[7:12], rna=rna)      #5-mer

	feat = pandas.DataFrame(featarray, index=data.index, columns=["Tm global_%s" % rna, "5mer_end_%s" %rna, "8mer_middle_%s" %rna, "5mer_start_%s" %rna])

	return feat

def gc_features(data):
	gc_count = data['30mer'].apply(countGC)
	gc_count.name = 'GC count'
	gc_above_10 = (gc_count > 10)*1
	gc_above_10.name = 'GC > 10'
	gc_below_10 = (gc_count < 10)*1
	gc_below_10.name = 'GC < 10'
	return gc_above_10, gc_below_10, gc_count

def normalize_features(data,axis):
	'''
	input: Pandas.DataFrame of dtype=np.float64 array, of dimensions
	mean-center, and unit variance each feature
	'''
	data -= data.mean(axis)
	data /= data.std(axis)
	# remove rows with NaNs
	data = data.dropna(1)
	if np.any(np.isnan(data.values)): raise Exception("found NaN in normalized features")
	return data

def apply_nucleotide_features(seq_data_frame, order, num_proc, include_pos_independent, max_index_to_use, prefix=""):

	fast = True

	if fast:
		feat_pd = seq_data_frame.apply(nucleotide_features, args=(order, include_pos_independent, max_index_to_use, prefix, 'pos_dependent'))
		feat_pi = seq_data_frame.apply(nucleotide_features, args=(order, include_pos_independent, max_index_to_use, prefix, 'pos_independent'))
	else:
		# old, much slower code
		feat_pd = pandas.DataFrame()
		feat_pi = pandas.DataFrame()
		for s in seq_data_frame.values:
			assert not (s==''), "string is empty"
			feat_pd_tmp, feat_pi_tmp = nucleotide_features(s, order, include_pos_independent, max_index_to_use, prefix=prefix)
			feat_pd = pandas.concat([feat_pd,feat_pd_tmp], axis=1)
			feat_pi = pandas.concat([feat_pi,feat_pi_tmp], axis=1)

		feat_pd = feat_pd.T
		feat_pi = feat_pi.T

	return feat_pd, feat_pi

def nucleotide_features(s, order, include_pos_independent, max_index_to_use, prefix="", feature_type='all'):
	'''
	compute position-specific order-mer features for the 4-letter alphabet
	(e.g. for a sequence of length 30, there are 30*4 single nucleotide features
		  and (30-1)*4^2=464 double nucleotide features
	'''
	if max_index_to_use is not None:
		s = s[:max_index_to_use]
	#assert(len(s)==30, "length not 30")
	#s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in, and were instructed to ignore in this way
	raw_alphabet = ['A', 'T', 'C', 'G']
	alphabet = ["".join(i) for i in itertools.product(raw_alphabet, repeat=order)]
	features_pos_dependent = np.zeros(len(alphabet)*(len(s)-(order-1)))
	features_pos_independent = np.zeros(np.power(len(raw_alphabet),order))


	#for position in range(0, len(s)-order, 1): JENN 9/4/2014 failing when len(s)=2
	for position in range(0, len(s)-order+1, 1):
		nucl = s[position:position+order]
		features_pos_dependent[alphabet.index(nucl)+(position*len(alphabet))] = 1.0
		features_pos_independent[alphabet.index(nucl)] += 1.0
	index_dependent = ['%s_pd.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_dependent))]

	if feature_type == 'all' or feature_type == 'pos_independent':
		if include_pos_independent:
			index_independent = ['%s_pi.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_independent))]
			if feature_type == 'all':
				return pandas.Series(features_pos_dependent,index=index_dependent), pandas.Series(features_pos_independent,index=index_independent)
			else:
				return pandas.Series(features_pos_independent, index=index_independent)

	if np.any(np.isnan(features_pos_dependent)): raise Exception("found nan features in features_pos_dependent")
	if np.any(np.isnan(features_pos_independent)): raise Exception("found nan features in features_pos_independent")

	return pandas.Series(features_pos_dependent, index=index_dependent)

def normalize_feature_sets(feature_sets):
	'''
	zero-mean, unit-variance each feature within each set
	'''

	print "Normalizing features..."
	t1 = time.time()

	new_feature_sets = {}
	for set in feature_sets:
		 new_feature_sets[set] = normalize_features(feature_sets[set],axis=0)
		 if np.any(np.isnan(new_feature_sets[set].values)):
			 raise Exception("found Nan feature values in set=%s" % set)

	t2 = time.time()
	print "\t\tElapsed time for normalizing features is %.2f seconds" % (t2-t1)

	return new_feature_sets

#############################
#                           #
#   Auxillary Functions     #
#                           #
#############################

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--indir','-i',help='filepath to directory which hosts BAM files',required=True)
	parser.add_argument('--filename','-n',help='name of BAM file',required=True)
	parser.add_argument('--fasta','-f',help='full filepath, including filename, of organism fasta file',required=True)
	parser.add_argument('--amino','-a',default=-1,help='amino acid information to be used in cutting efficiency.'
													   ' Default -1 indicating to not use amino acid information')
	parser.add_argument('--peptide','-p',default=-1,help='peptide information to be used in cutting efficiency. '
														 'Default -1 indicating to not use peptide information')
	parser.add_argument('--nopos','-q',default='../Rule_Set_2_scoring/saved_models/V3_model_nopos.pickle',
						help='filepath to V3_model_nopos.pickle. Defaults to model saved in package directory')
	parser.add_argument('--full','-w',default='../Rule_Set_2_scoring/saved_models/V3_model_full.pickle',
						help='filepath to V3_model_full.pickle. Defaults to model saved in package directory')

	args = parser.parse_args()
	indir = args.indir
	filename = args.filename
	fasta = args.fasta
	aa_cut = args.amino
	per_peptide = args.peptide
	V3_model_nopos = args.nopos
	V3_model_full = args.full

	return indir,filename,fasta,aa_cut,per_peptide,V3_model_nopos,V3_model_full

def serial_on_target_doench(seq,aa_cut,per_peptide,model):
	if seq[25:27] == 'GG':
		return predict(seq, aa_cut, per_peptide, model=model)
	else:
		sys.stderr.write('WARNING: GG PAM not encountered')

def donech_scoring_parameters(aa_cut,per_peptide,V3_model_nopos,V3_model_full):
	"""Load the correct Donech scoring model

	Convention: The V3_model_nopos.pickle is the model which does not utilize amino acid or peptide information
	in predicting a guideRNA's cutting efficiency. To access this model aa_cut and per_petide are both set to -1.
	However, this model will be called if either aa_cut or per_peptide is set to -1. The V3_model_full.pickle is
	the model which utilizes amino acid and peptide information in its prediction of guideRNA cutting efficiency.
	This model will be call if both aa_cut and per_peptide are not set to -1. V3_model_nopos and V3_model_full
	are the full filepaths to these models.

	Args:
	aa_cut: indicates whether amino acid information should be used in predicting guideRNA cutting efficiency. If
			set to -1 then V3_model_nopos.pickle will be loaded

	per_peptide: indicates whether peptide information should be used in predicting guideRNA cutting efficiency. If
				 set to -1 then The V3_model_full.pickle will be loaded.

	V3_model_nopos: full filepath to V3_model_nopos.pickle

	V3_model_full: full filepath to  V3_model_full.pickle
	"""
	model_file_1 = V3_model_nopos
	model_file_2 = V3_model_full
	if (aa_cut == -1) or (per_peptide == -1):
	   model_file = model_file_1
	else:
		model_file = model_file_2
	try:
		with open(model_file, 'rb') as f:
			model = pickle.load(f)
	except:
		raise Exception("could not find model stored to file %s" % model_file)
	return model

def bamfile_info_extraction(indir,filename):
	"""Extract chromosome, coordinate, strand, original sequence from BAM file

	Args:
	indir = directory to where BAM file is hosted. This is filepath does not
			include the name of the file itself.
	filename = name of the BAM file from which information will be extracted

	Note: the functionality comes from a system call to Samtools and the AWK
		  language.
	"""
	infile = indir + '/' + filename
	outfile = indir + '/' + 'CutEffScoreExtract_' + filename.replace('.bam','.txt')

	#extract relevant information from BAM file
	cmd1 = "samtools view " + infile + " | awk '{print $3,$4,$2,$1}' - > " + outfile
	os.system(cmd1)

	return ('cutting efficiency data for %s extracted' % (filename)),outfile

def fasta(fasta_file):
	"""Load organism fasta file for use in pyfaidx module

	Args:
	fasta_file = the full filepath, including the file itself, to the organism's fasta file

	Note: an index of the fasta file should also be present in the same directory. This can
		  be produced using samtools faidx command and will have the suffix .fai
	"""
	org = Fasta(fasta_file)
	return ('%s accessed' % fasta_file),org

def cutting_efficiency_w_doench_rule2(outfile,indir,filename,org,aa_cut,per_peptide,model):
	"""Compute cutting efficiency using Rule 2 model

	Args:
	outfile = output file of bamfile_info_extraction()
	indir = filepath to directory hosting BAM file. Same input as given to bamfile_info_extraction()
	filename = name of BAM file in indir which will be assessed. Same input as given to bamfile_info_extraction()
	org = output of fasta()
	aa_cut = same input as given to donech_scoring_parameters()
	per_peptide = same input as given to donech_scoring_parameters()
	model = output of donech_scoring_parameters()
	"""
	bam_extracts = open(outfile,'r')
	cutting_efficiency_outfile = open(indir + '/' + 'CutEffScore_' + filename.replace('.bam','.txt'),'w')
	for line in bam_extracts:
		clean_line = line.lstrip().rstrip()
		parts = clean_line.split()
		chrom,start_coord,strand,sequence = parts[0],int(parts[1]),int(parts[2]),parts[3]
		if strand == 0:
			start_coord_30 = start_coord - 5
			end_coord_30 = start_coord + len(sequence) + 2
			try:
				r20_to_30mer = str(org[chrom][start_coord_30:end_coord_30]).upper()

				if len(r20_to_30mer) == 30:

					if 'N' not in r20_to_30mer:
						CutEffScore = 'ds:Z:' + str(serial_on_target_doench(r20_to_30mer,aa_cut,per_peptide,model))
						cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																								  end_coord_30,strand,CutEffScore))
					else:
						CutEffScore = 'ds:Z:' + str(serial_on_target_doench('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',aa_cut,per_peptide,model))
						cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																								  end_coord_30,strand,CutEffScore))
				else:
					cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																							  end_coord_30,strand,'ds:Z:-1'))
			except:
				cutting_efficiency_outfile.write(
					'%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence, chrom, start_coord_30,
															 end_coord_30, strand, 'ds:Z:-1'))

		if strand == 16:
			start_coord_30 = start_coord - 4
			end_coord_30 = start_coord + len(sequence) + 3
			try:
				r20_to_30mer = str(Seq.Seq(str(org[chrom][start_coord_30:end_coord_30]).upper()).reverse_complement())
				if len(r20_to_30mer) == 30:

					if 'N' not in r20_to_30mer:
						CutEffScore = 'ds:Z:' + str(serial_on_target_doench(r20_to_30mer,aa_cut,per_peptide,model))
						cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																								  end_coord_30,strand,CutEffScore))
					else:
						CutEffScore = 'ds:Z:' + str(serial_on_target_doench('GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',aa_cut,per_peptide,model))
						cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																								  end_coord_30,strand,CutEffScore))
				else:
					cutting_efficiency_outfile.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence,chrom,start_coord_30,
																							  end_coord_30,strand,'ds:Z:-1'))
			except:
				cutting_efficiency_outfile.write(
					'%s \t %s \t %s \t %s \t %s \t %s \n' % (sequence, chrom, start_coord_30,
															 end_coord_30, strand, 'ds:Z:-1'))

	cutting_efficiency_outfile.close()
	bam_extracts.close()
	return ('cutting efficiency for %s computed' % (filename))

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

def v3_model_path_fidelity(V3_model_nopos,V3_model_full):
		exit_code = 0

		if str(V3_model_nopos).endswith('V3_model_nopos.pickle'):
			exit_code += 0
		else:
			exit_code += 1
			sys.stderr.write('ERROR: filepath does not end in V3_model_nopos.pickle')

		if str(V3_model_full).endswith('V3_model_full.pickle'):
			exit_code += 0
		else:
			exit_code += 1
			sys.stderr.write('ERROR: filepath does not end in V3_model_full.pickle')

		if exit_code == 0:
			return 0
		else:
			return 1

#####################
#                   #
#   Core Function   #
#                   #
#####################

def core(indir,filename,fasta_file,aa_cut,per_peptide,V3_model_nopos,V3_model_full):

	#indir,filename,fasta_file,aa_cut,per_peptide,V3_model_nopos,V3_model_full = arg_parser()

	#ensure V3 models end with model file
	v3_model_paths = v3_model_path_fidelity(V3_model_nopos,V3_model_full)

	#ensure file existence
	db_exist = file_path_verification(indir + '/' + filename)
	fasta_exist = file_path_verification(fasta_file)
	v3_nopos_exist = file_path_verification(V3_model_nopos)
	v3_full_exist = file_path_verification(V3_model_full)

	if db_exist or fasta_exist or v3_nopos_exist or v3_full_exist or v3_model_paths != 0:
		sys.exit()

	else:

		model = donech_scoring_parameters(aa_cut,per_peptide,V3_model_nopos,V3_model_full)

		msg,outfile = bamfile_info_extraction(indir,filename)
		sys.stdout.write(msg + '\n')

		msg,org = fasta(fasta_file)
		sys.stdout.write(msg + '\n')

		msg = cutting_efficiency_w_doench_rule2(outfile,indir,filename,org,aa_cut,per_peptide,model)
		sys.stdout.write(msg + '\n')

		cmd2 = 'rm ' + outfile
		os.system(cmd2)
		sys.stdout.write('clean up' + '\n')
		return

#if __name__ == '__main__':
#    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
