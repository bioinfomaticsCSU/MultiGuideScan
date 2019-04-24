__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""computes the CFD and Hsu et al. gRNA specificity score

Note: Cutting specificity scores are computed using the CFD method devised by Doench et al. Nature Biotechnology 2016.
Functions under --CFD Functions-- are used wholesale or with slight modification from original paper. Functions under
--Auxillary Functions-- are written by authors. Author written functions are summarized below:

#####################################
#                                   #
#   Auxillary Functions summaries   #
#                                   #
#####################################

file_path_verification(database) = function that verifies that the provided filepath exists before further processing occurs
hsu_scoring(string1, string2) = computes Hsu score for a gRNA
fasta(fasta_file) = Load organism fasta file for use in pyfaidx module
bamfile_info_extraction(indir, filename) = Extract chromosome, coordinate, strand, original sequence from BAM file
specificity_score(outfile,bamfile,org,mm_scores,pam_scores,scoring_outfile) = computes CFD and/or Hsu specificity score for gRNA
organism_chromosome_identity_and_length_confirmation(bam_file) = creates dictionary object containing chromosome numbers and lengths
query_bam(guidesbam, coords, offcoords=True, onebased=False) = query BAM file containing guideRNA database.
get_bed_lines(guides, offtargets=True) = Create lines in BED format form a list of GuideRNA objects

"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import sys
import ast
import timeit
import pickle
import numpy as np

from Bio import Seq
from pyfaidx import Fasta
import pysam

import util

#####################
#                   #
#   CFD Function    #
#                   #
#####################

def calc_cfd(wt,sg,pam,mm_scores,pam_scores):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            try:
                key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
                score*= mm_scores[key]
            except KeyError:
                continue
    score*=pam_scores[pam]
    return (score)

def get_mm_pam_scores(mms,pams):
    try:
        mm_scores = pickle.load(open(mms,'rb'))
        pam_scores = pickle.load(open(pams,'rb'))
        return (mm_scores,pam_scores)
    except:
        raise Exception("Could not find file with mismatch scores or PAM scores")

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#############
#           #
#   Class   #
#           #
#############

class GuideRNA:

    """CRISPR-Cas system guideRNA with additional info including off-targets"""

    def __init__(self, seq, coord, score=-1, offdist=-1, maxoffcount=-1,
                 offtargets=[]):
        """Initialize guideRNA with all relevant info.

        Args:
        seq: str with guideRNA sequence including PAM
        coord: tuple (chromosome, start, end, strand)
               where chromosome is str, start and end are int,
               strand is either '+' or '-';
               start and end are positions of first and last characters
               of guideRNA in the genome (i.e., both stard and end
               are part of guideRNA);
               assume coordinates are 0-based
        score: cutting efficiency score for the guideRNA; float between 0
               (low cutting efficiency) and 1 (high cutting efficiency);
               if -1, score is undefined
        offdist: int, shows what offdist parameter was used
                 when calculating off-target info;
                 values larger than offdist cannot appear
                 as Hamming distance in off-target info;
                 in particular, empty off-target info list
                 can have different meaning for different offdist values;
                 if -1, no assumptions can be made about offtargets
        maxoffcount: int, shows what maxoffcount parameter was used
                     when calculating off-target info;
                     if -1, no assumptions can be made about offtargets
        offtargets: list of tuples corresponding to off-target k-mers
                    (Hamming distance (int) from the k-mer to the guideRNA,
                     number of occurrences (int) of this k-mer in the genome,
                     list of genomic coordinates where this k-mer occurs);
                    the list of coordinates, if not empty, consists of tuples
                    (chrom (str),
                     0-based coord of first position of occurrence (int),
                     strand (str))
        """
        self.seq = seq
        self.coord = coord
        self.score = score
        self.offdist = offdist
        self.maxoffcount = maxoffcount
        self.offtargets = offtargets

    def __repr__(self):
        """String representation. Currently same as self.bed()."""
#         osummary = self.offtarget_summary()
#         offtarget_summary_str = '|'.join('%s:%s' % (i, count)
#                                          for i,count in enumerate(osummary)
#                                          if i >= 2)
# #        score = self.offtarget_score(osummary)
#         s = '%s:%s-%s:%s\t%s\t%s\t%s' % \
#             (self.coord + (self.seq, sum(osummary), offtarget_summary_str))
#         return s
        return self.bed()

    def bed(self, detailed=False):
        """String representation in (extended) BED format.

        Args:
        detailed: if True, add additional fields; useful if the
                  resulting line appears in the same file
                  with off-targets and need to distinguish
                  guideRNAs from off-targets
        """
        fields = self.bed_fields()
        s = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % fields[:9]
        if detailed:
            s += '\t%s\t%s' % fields[9:]
        return s

    def bed_fields(self):
        """Produce tuple of fields to use in (extended) BED format.

        Current convention for BED format representation is stored
        as variable BED_FORMAT_DESCRIPTION.

        BED format is defined here:
        https://genome.ucsc.edu/FAQ/FAQformat.html#format1
        According to BED format specification, output coordinates are 0-based
        and end coordinate is not part of guideRNA.
        """
        osummary = self.offtarget_summary()
        offtarget_summary_str = '|'.join('%s:%s' % (i, count)
                                         for i,count in enumerate(osummary)
                                         if i >= 2)
        offtargscore = self.offtarget_score(osummary)
        fields = (self.coord[0], self.coord[1], self.coord[2] + 1, self.seq,
                  int(100 * self.score) if 0 <= self.score <= 1 else '*',
                  self.coord[3], sum(osummary),
                  offtarget_summary_str,
                  int(100 * offtargscore) if 0 <= offtargscore <= 1 else '*',
                  0, '*')
        return fields

    def offtarget_bed_fields(self):
        """Produce list of BED tuples for all off-targets with known coords.

        Current convention for BED format representation is stored
        as variable BED_FORMAT_DESCRIPTION.
        """
        l = []
        for hammdist, count, offcoords in self.offtargets:
            for chrom, coord, strand in offcoords:
                if strand == '+':
                    start, end = coord, coord + len(self.seq)
                elif strand == '-':
                    start, end = coord - len(self.seq) + 1, coord + 1
                fields = (chrom, start, end, self.seq, '*', strand,
                          '*', '*', '*', hammdist, '*')
                l.append(fields)
        return l

    def offtarget_summary(self, maxdist=None):
        """Find off-target count for each Hamming distance.

        Args:
        maxdist: maximum Hamming distance to take into account;
                 by default, use all values for which off-target
                 info is known (i.e., up to offdist)

        Return:
        tuple (count of off-targets for Hamming distance 0,
               count of off-targets for Hamming distance 1,
               count of off-targets for Hamming distance 2, etc.)
        """
        maxdist = maxdist or self.offdist
        if self.offdist != -1:
            return tuple(sum(p[1] for p in self.offtargets if p[0] == i)
                         for i in xrange(maxdist + 1))

    def offtarget_score(self, osummary=None):
        """Score off-target potential.

        Assign score from 0 to 1. 1 means little off-target potential, 0 means
        a lot of off-targets expected.
        Score currently is very naive:
        1 / (1 + total weighted sum of off-targets)
        where each off-target at Hamming distance D contributes to
        the weighted sum with weight 1 / D

        Args:
        osummary: output of offtarget_summary(); if None, compute from scratch
        """
        if not osummary:
            osummary = self.offtarget_summary()
        score = 1.0 / (1 + sum(v / float(i) for i,v in enumerate(osummary)
                               if i > 0))
        return score

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def kmer_exact_occurrence_dictionary(kmer_counts_file):
    """creates dictionary object which stores kmers as key and kmer occurrence in genome as values

    Input:
    kmer_counts_file: absolute filepath to XXX_all_kmers_counted.txt file

    """
    kmer_dictionary = {}
    with open(kmer_counts_file, 'r') as kmers:
        for line in kmers:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split()
            if kmer_dictionary.has_key(parts[1]):
                sys.stdout.write('kmer duplication detected %s %s \n' % (parts[0], parts[1]))
            else:
                kmer_dictionary[parts[1]] = parts[0]

        sys.stdout.write('kmer dictionary generated \n')

    return kmer_dictionary

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

def fasta(fasta_file):
	"""Load organism fasta file for use in pyfaidx module

	Args:
	fasta_file = the full filepath, including the file itself, to the organism's fasta file

	Note: an index of the fasta file should also be present in the same directory. This can
		  be produced using samtools faidx command and will have the suffix .fai
	"""
	org = Fasta(fasta_file)
	return ('%s accessed' % fasta_file),org

def bamfile_info_extraction(indir, filename):
    """Extract chromosome, coordinate, strand, original sequence from BAM file

	Args:
	indir = directory to where BAM file is hosted. This is filepath does not
			include the name of the file itself.
	filename = name of the BAM file from which information will be extracted

	Note: the functionality comes from a system call to samtools and the AWK
		  language.
	"""
    infile = indir + '/' + filename
    outfile = indir + '/' + 'CutEffScoreExtract_' + filename.replace('.bam', '.txt')

    # extract relevant information from BAM file
    cmd1 = "samtools view " + infile + " | awk '{print $3,$4,$2,$1}' - > " + outfile
    os.system(cmd1)

    return ('cutting efficiency data for %s extracted' % (filename)), outfile

def specificity_score(outfile, bamfile, org, mm_scores, pam_scores, scoring_outfile,kmer_dictionary):
    """computes CFD and/or Hsu specificity score for gRNA

    Args:
    outfile: output file of bamfile_info_extraction()
    bamfile: absolute filepath to database BAM file
    org: second output object of fasta(fasta_file)
    mm_scores: path to CFD mismatch score pickle
    pam_scores: path to CFD pam score pickle
    scoring_outfile: absolute filepath to output file where specificity scores will be written

    """
    gRNA_dictionary = {}
    bam_extracts = open(outfile, 'r')
    scoring_file = open(scoring_outfile, 'w')
    # scoring_file.write('%s \t %s \t %s \n' % ('unique_identifier','mean_cfd_score','sum cfd score'))
    count = 0
    for line in bam_extracts:
        clean_line = line.lstrip().rstrip()
        parts = clean_line.split()
        count += 1
        chrom, start_coord, strand_db, sequence_db = parts[0], int(parts[1]), int(parts[2]), parts[3]
        if strand_db == 16 or strand_db == 0:

            if strand_db == 0:
                start_coord_30 = start_coord
                end_coord_30 = start_coord + len(sequence_db)
                read = '%s:%s-%s' % (chrom, start_coord_30 - 1, end_coord_30)

            elif strand_db == 16:
                start_coord_30 = start_coord
                end_coord_30 = start_coord + len(sequence_db)
                read = '%s:%s-%s' % (chrom, start_coord_30-1, end_coord_30-1) #maybe shift end coord


            guides = query_bam(bamfile,read,offcoords=True)
            if guides:
                for j in range(len(guides)):
                    sequence = guides[j].seq
                    start_coord_30,end_coord_30 = guides[j].bed().split()[1],guides[j].bed().split()[2]
                    if gRNA_dictionary.has_key(sequence):
                        continue
                    else:
                        gRNA_dictionary[sequence] = guides[j].offtarget_score()
                        strand_wo =  guides[j].coord[3]
                        offtargets = guides[j].offtargets
                        offtarget_lst = []
                        if offtargets:
                            for i in offtargets:
                                #print i
                                chromosome = i[2][0][0]
                                strand = i[2][0][2]
                                if strand == '+':
                                    start_coord = i[2][0][1]
                                    end_coord = start_coord + len(sequence)
                                    ot_sequence = list(str(org[chromosome][start_coord:end_coord]).upper())
                                    ot_sequence[-3] = 'N'
                                    ot_sequence = ''.join(ot_sequence)

                                elif strand == '-':
                                    end_coord = i[2][0][1]
                                    start_coord = end_coord - len(sequence)
                                    ot_sequence = list(str(
                                        Seq.Seq(
                                            str(org[chromosome][start_coord + 1:end_coord + 1]).upper()).reverse_complement()))
                                    ot_sequence[-3] = 'N'
                                    ot_sequence = ''.join(ot_sequence)

                                pam = ot_sequence[-2:]
                                sg = ot_sequence[:-3]
                                ot_sequence_occurence = kmer_dictionary[ot_sequence] #enumerate occurrence of OT
                                cfd_score = calc_cfd(sequence, sg, pam, mm_scores, pam_scores)
                                total_cfd_contribution = cfd_score * float(ot_sequence_occurence) #CFD w occurrence
                                offtarget_lst.append(total_cfd_contribution)

                                #print ot_sequence, sequence, cfd_score
                                count, mismatches = 0, 0
                                for chr1, chr2 in zip(ot_sequence, sequence):
                                    if chr1 != chr2:
                                        mismatches += 1
                                        count += 1
                                    else:
                                        count += 1

                                #print ot_sequence, sequence, strand, i[0], mismatches,cfd_score

                            offtarget_array = np.array(offtarget_lst)
                            cfd_aggregate_score = 1.0 / (1.0 + offtarget_array.sum())
                            #print 'total aggretate score for ', sequence, cfd_aggregate_score

                            #print ot_sequence,sequence,strand,cfd_aggregate_score,'tag'
                            scoring_file.write('%s \t %s \t %s \t %s \t %s \t cs:Z:%.8f \n' % (
                                sequence, chrom, start_coord_30, end_coord_30, strand_wo, cfd_aggregate_score))  # ot_np.sum()
                        else:
                            #print 'passed, no ot',sequence
                            scoring_file.write('%s \t %s \t %s \t %s \t %s \t cs:Z:%.8f \n' % (
                                sequence, chrom, start_coord_30, end_coord_30, strand_wo, 1))
            else:
                continue
        else:
            continue

    scoring_file.close()
    bam_extracts.close()
    print count
    return ('cutting efficiency for %s computed' % (outfile))

def organism_chromosome_identity_and_length_confirmation(bam_file):
    """creates dictionary object containing chromosome numbers and lengths

    Input:
    bam_file: absolute filepath to BAM database

    """
    verification_dict = {}
    bamfile = pysam.AlignmentFile(bam_file,'rb')
    chrom_data = bamfile.header
    for i in range(len(chrom_data['SQ'])):
        verification_dict[chrom_data['SQ'][i]['SN']] = chrom_data['SQ'][i]['LN']
    bamfile.close()
    return verification_dict

def query_bam(guidesbam, coords, offcoords=True, onebased=False):
    """Query BAM file containing guideRNA database.

    Input:
    guidesbam: path to BAM file
    coords: genomic coordinates, e.g., str 'chrX:3364088-3372035'
            or tuple (chr, start, end), e.g., ('chrX', 3364088, 3372035) ;
            assume coordinates are 0-based
    offcoords: if True, parse off-target coordinates (takes time), otherwise
               keep them in int encoded form (can be decoded later
               with util.map_int_to_coord());
               if parsed, coordinates of off-targets are 0-based
    onebased: if True, input coordinates are 1-based (default is 0-based)

    Output:
    list of instances of class GuideRNA
    """

    #create dictionary to verify query coordinates correspond to organism chromosome identity and length
    #chromosome_dictionary = organism_chromosome_identity_and_length_confirmation(guidesbam)

    if isinstance(coords, basestring):
        coords = coords.strip().split(':')
        coords = (coords[0],) + tuple(int(x) for x in coords[1].split('-'))
    chrom, start, end = coords

    #confirm coordinate exists in genome
    """
    if chromosome_dictionary.has_key(chrom):
        if start <= chromosome_dictionary[chrom] and end <= chromosome_dictionary[chrom]:
            pass
        else:
            sys.stderr.write('WARNING: %s and/or %s fall off %s which has length %s \n' % (start,end,chrom,chromosome_dictionary[chrom]))
            return
    else:
        sys.stderr.write('WARNING: %s is not a recognized chromosome for this organism \n' % (chrom))
        return
    """
    if onebased:
        start -= 1 # switch to 0-based coordinates, as expected by pysam
        end -= 1   # and used internally in our pipeline (in the trie etc.)
    bam = pysam.AlignmentFile(guidesbam, "rb")
    guides = []
    argsdict = ast.literal_eval(bam.header['CO'][3])
    genome = argsdict['genome']
    delim = util.get_nonexist_int_coord(genome)
    for r in bam.fetch(chrom, start, end):
        seq = r.qname
        positions = r.get_reference_positions()
        # note: pysam returns 0-based coordinates
        r_start, r_end = min(positions), max(positions)
        if r_start < start or r_end > end:
            # ignore guideRNAs not completely within query genomic region
            continue
        strand = '+' if not r.is_reverse else '-'
        tags = dict(r.get_tags())
        offtargets = tags.get('of', [])
        offdist = tags.get('od', -1)
        # if 'ds' in r.tags[len(r.tags)-1]:
        #     score = r.tags[len(r.tags)-1][1]
        score = tags.get('ds', -1)

        try:
            if util.is_number(score):
                score = float(score)
            else:
                score = -1
        except ValueError:
            sys.stderr.write('ERROR: score %s could not be converted to float, skipping region' % (score))
            return

        if offdist == -1:
            offdist = argsdict.get('offdist', -1)
        maxoffcount = tags.get('oc', -1)
        if maxoffcount == -1:
            maxoffcount = argsdict.get('maxoffcount', -1)
        if isinstance(offtargets, basestring):
            offtargets = util.hex_to_offtargetinfo(offtargets, delim,
                                                   offcoords)
            offtargets = [(p[0],
                           (1 if p[1] == 0 else p[1]) if not offcoords
                           else (1 if p[1][0] == 0 else p[1][0]),
                           [] if not offcoords
                           else [util.map_int_to_coord(x, genome,
                                                       onebased=False,
                                                       strout=False)
                                 for x in p[1][1:]])
                          for p in offtargets]
        guiderna = GuideRNA(seq=seq, coord=(chrom, r_start, r_end,
                                            strand),
                            score=score,
                            offdist=offdist, maxoffcount=maxoffcount,
                            offtargets=offtargets)
        guides.append(guiderna)
    bam.close()
    return guides

def get_bed_lines(guides, offtargets=True):
    """Create lines in BED format form a list of GuideRNA objects.

    Args:
    guides: list of GuideRNA objects such as query_bam() output
    offtargets: if True, also include a line for each known off-target
                of each guideRNA from guides
    """
    bedlines = []
    for guide in guides:
        bedlines.append(guide.bed_fields())
        if offtargets:
            bedlines.extend(guide.offtarget_bed_fields())
    return bedlines

#####################
#                   #
#   Core Function   #
#                   #
#####################

def core(indir,filename,fasta_file,mismatch_score,pam_score,kmer_counts_file):

    """
    indir = '/Users/pereza1/Projects/Ventura/CRISPR/data/mm10/test'
    filename = 'test.bam'
    bamfile = '%s/%s' % (indir, filename)
    fasta_file = '/Users/pereza1/Reference_Genomes/mm10/mm10.fa'
    mismatch_score = '%s/mismatch_score.pkl' % ('/Users/pereza1/Projects/Ventura/CRISPR/code/auxillary/CFD_Scoring')
    pam_score = '%s/pam_scores.pkl' % ('/Users/pereza1/Projects/Ventura/CRISPR/code/auxillary/CFD_Scoring')
    """

    #database file
    bamfile = '%s/%s' % (indir, filename)

    #cfd scores
    mm_scores,pam_scores = get_mm_pam_scores(mismatch_score,pam_score)

    #fasta index
    msg,org = fasta(fasta_file)
    sys.stdout.write(msg + '\n')

    #bam file data extraction
    msg,outfile = bamfile_info_extraction(indir, filename)
    sys.stdout.write(msg + '\n')

    #scoring outfile
    scoring_outfile = '%s/scoring_outfile_%s' % (indir,filename.replace('bam','txt'))

    #kmer counts dictionary
    kmer_dictionary = kmer_exact_occurrence_dictionary(kmer_counts_file)

    #scoring
    start = timeit.default_timer()
    msg = specificity_score(outfile,bamfile,org,mm_scores,pam_scores,scoring_outfile,kmer_dictionary)
    sys.stdout.write(msg + '\n')
    end = timeit.default_timer()
    sys.stdout.write('time to complete is %s seconds\n' % (end-start))
    return
