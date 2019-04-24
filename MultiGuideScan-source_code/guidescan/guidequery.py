__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""User interface to precomputed guideRNA library stored in BAM format

GuideRNA: class to store and show guideRNAs with all relevant info
query_bam(): query BAM file containing guideRNA database
"""

#########################
#                       #
#   Set Descriptions    #
#                       #
#########################

# description of BED format for guideRNAs or their off-targets
BED_FORMAT_DESCRIPTION = """
        0: chromosome name;
        1: start coordinate on chromosome (0-based);
        2: end coordinate on chromosome (not part of guideRNA or off-target);
        3: sequence of guideRNA + PAM (may include N for unknown nucleotide);
        4: cutting efficiency score (between 0 and 100, 100 is best
           efficiency) for guideRNA, "*" for off-targets or if unknown;
        5: cutting specificity score for gRNA, * if unknown
        6: strand;
        7: total number of off-targets of guideRNA, undefined for off-targets;
        8: summary of off-targets for a guideRNA in the form
           "#mismatches:#offtargets" separated by "|" ("*" if unknown);
        9: guideRNA off-target score (between 0 and 100, 100 means no
           off-targets), "*" for off-targets or if unknown;
        10: number of mismatches of genomic sequence with the guideRNA+PAM
           (always 0 for guideRNAs, some value larger than 0 for off-targets);
        11: sequence in the genome ("*" if unknown), for guideRNA will coincide
            with guideRNA+PAM sequence (field 3) with all N substituted
            to actual sequence content, for off-targets will be distant from
            this guideRNA+PAM sequence by a certain number of mismatches
            (field 8);
"""
BED_HEADER = ('chrom', 'chromStart', 'chromEnd', 'guideRNA+PAM',
              'cutEffScore', 'strand',
              'offtargCount', 'offtargSummary', 'offtargScore',
              'numMismatch', 'seq')

BED_HEADER_BED_INPUT_WITHIN = ('chromosome', 'start', 'end', 'target_sequence', 'cutting_efficiency_score','cutting_specificity_score','strand',
                               'offtargets_sum', 'offtarget_summary', 'quality_score', 'annotation',
                               'gRNA label')

BED_HEADER_BED_INPUT_FLANKS = ('chromosome', 'start', 'end', 'target_sequence',
                               'cutting_efficiency_score','cutting_specificity_score', 'strand', 'offtargets_sum', 'offtarget_summary',
                               'quality_score',
                               'x', 'left_flank_distance', 'annotation', 'left_gRNA_label', '', '', 'chromosome',
                               'start', 'end', 'target_sequence', 'cutting_efficiency_score','cutting_specificity_score',
                               'strand', 'offtargets_sum', 'offtarget_summary', 'quality_score', 'x',
                               'right_flank_distance',
                               'annotation', 'right_gRNA_label', '', 'total_flanking_distance', 'total_deletion_size')

#################
#               #
#   Libraries   #
#               #
#################

import re
import os
import sys
import ast
import csv
import argparse
import subprocess
import xlwt

import pysam
from itertools import repeat
from collections import defaultdict

import util

#############
#           #
#   Class   #
#           #
#############

class GuideRNA:

    """CRISPR-Cas system guideRNA with additional info including off-targets"""

    def __init__(self, seq, coord, score=-1, offdist=-1, maxoffcount=-1, specificity=-1,
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
        specificity: cutting specificity score, how likely gRNA is to hit on-target, given information on the gRNA's
                        offtargets
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
        self.specificity = specificity

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
                  int(100 * self.score) if 0 <= self.score <= 1 else '*',self.specificity if 0 <= self.specificity else '*',
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

def ucsc_compatibility_insurance(guides):
    """make coordinates compatible with UCSC coordinates

    Args:
        guides: GuideRNA object which contains sgRNAs. Output of guidequery.query_bam

    Returns: GuideRNA object with coordinates made compatible with UCSC genomic coordinates

    """
    try:
        bedlines = get_bed_lines(guides, offtargets=False)

        # capture data needed to reconstruct GuideRNA class
        offdist = guides[0].offdist
        maxoffcount = guides[0].maxoffcount

        # off-target information captured for each sgRNA before class break
        offtarget_dictionary_GuideRNA = defaultdict(list)
        for j in range(len(guides)):
            key, value = str(guides[j]).split()[3], guides[j].offtargets
            offtarget_dictionary_GuideRNA[key].append(value)

        # reconstruct GuideRNA class
        bedlines_lst = []
        for i in range(len(bedlines)):
            if type(bedlines[i][4]) != int:
                if bedlines[i][6] == '+':
                    start = bedlines[i][1] + 1
                    end = bedlines[i][2] - 2
                elif bedlines[i][6] == '-':
                    start = bedlines[i][1] + 1
                    end = bedlines[i][2] - 2

                seq, coord, score, specificity = bedlines[i][3], (
                    bedlines[i][0], start, end, bedlines[i][6]), -1, bedlines[i][5]
            else:
                if bedlines[i][6] == '+':
                    start = bedlines[i][1] + 1
                    end = bedlines[i][2] - 2
                elif bedlines[i][6] == '-':
                    start = bedlines[i][1] + 1
                    end = bedlines[i][2] - 2

                seq, coord, score, specificity = bedlines[i][3], (
                    bedlines[i][0], start, end, bedlines[i][6]), \
                                                 (bedlines[i][4] / 100.0), bedlines[i][5]

            for q in offtarget_dictionary_GuideRNA[seq]:
                select_guide = GuideRNA(seq=seq, coord=coord, score=score, specificity=specificity, offdist=offdist,
                                        maxoffcount=maxoffcount,
                                        offtargets=q)
                bedlines_lst.append(select_guide)

        return bedlines_lst

    except IndexError:
        sys.stderr.write('WARNING: Index Error detected, returning original GuideRNA object without ucsc coordinate modification\n')
        return guides

def file_existance_fasta_index(input_file):
    """Ensures a file exists at the stated path

    Input:
    input_file: absolute filepath of the file to be processed

    """
    if os.path.exists('%s.fai' % input_file) and os.path.exists(input_file):
        return 0
    else:
        return 1

def blat_processing(blat, fasta, infile, outdir, verbose=True):
    """

    :param blat: absolute filepath to blat binary, (unix command line: which blat)
    :param fasta: absolute filepath to indexed fasta file
    :param infile: input file with sequences to process; fasta format
    :param outdir: absolute filepath to output directory
    :param verbose: return stdout
    :return: returns absolute filepath of blat coordinates file

    """
    blat_check = blat_local_installation_verification(blat)

    if blat_check == 0:
        out = '%s/out_for_blat.txt' % (outdir)
        run_blat(blat, fasta, infile, out)
        outfile = blat_coordinate_file(out, verbose)
        os.remove(out)
        return outfile
    else:
        sys.stderr.write('blat not run on %s due to blat not being locally installed on system\n' % infile)
        return

def blat_coordinate_file(out,verbose=False):
    """generates a tab-delimited file with genomic coordinates for GuideScan

    Input:
    out: first output object of run_blat()
    verbose: write out lines in blat file

    Output:
    a tab-delimited file with the best blat estimate of the genomic coordinates for the input sequences

    Note:
    this function will only return coordinates for a sequence, if that sequence has a perfect alignment, according to BLAT,
    in a given genome.

    """

    outdir = '/'.join(out.split('/')[:-1])
    with open('%s/blat_coordinates.txt' % (outdir),'w') as outfile:
        with open(out,'r') as infile:
            for line in infile:
                clean_line = line.lstrip().rstrip()
                parts = clean_line.split()
                if verbose == True:
                    print clean_line
                else:
                    pass
                if parts:
                    try:
                        score,qsize = int(parts[0]),int(parts[10])
                        if score == qsize:
                            blat_coordinate = '%s:%s-%s' % (parts[13], parts[15], parts[16])
                            outfile.write('%s\n' % blat_coordinate)
                    except ValueError:
                        continue

            sys.stdout.write('coordinates file written to %s/blat_coordinates.txt\n' % outdir)

    return '%s/blat_coordinates.txt' % (outdir)

def run_blat(blat,fasta,infile,outfile):
    """executes blat command

    Input:
    blat: absolute filepath of blat
    fasta: absolute filepath to indexed fasta file for organism under consideration
    infile: input fasta file
    outfile: output fasta file

    Output:
    blat_file: output file from blat

    """
    cmd = '%s %s %s %s' % (blat,fasta,infile,outfile)
    os.system(cmd)
    sys.stdout.write('blat run complete\noutput written to %s\n' % (outfile))

def blat_local_installation_verification(blat):
    """verifies that blat is locally installed on system

    Inputs:
    blat: absolute filepath of blat

    Outputs:
    returns 0 if blat is locally installed and 1 if not

    """
    ucsc_executables = 'http://hgdownload.soe.ucsc.edu/admin/exe/'
    try:
        check = os.system(blat)
        if check == 65280:
            sys.stdout.write('blat is locally installed\n')
            return 0
        else:
            sys.stderr.write('blat is not locally installed at %s, please install from %s\n' % (blat, ucsc_executables))
            return 1
    except NameError:
        sys.stderr.write('blat is not locally installed at %s, please install from %s\n' % (blat,ucsc_executables))
        return 1
        #sys.exit(1)

def gene_dictionary_from_bed_file(bed_file):
    """generates a gene dictionary object where gene name is key and coordinates are values

    Input:
    bed_file: absolute filepath to a bed file with gene annotations

    Output:
    gene_dictionary: dictionary object with gene name as key and genomic coordinate as value

    """
    gene_dictionary = {}
    with open(bed_file,'r') as bed:
        for line in bed:
            clean_line = line.lstrip().rstrip()
            parts = clean_line.split()
            gene_dictionary[parts[3]] = '%s:%s-%s' % (parts[0],parts[1],parts[2])

    sys.stdout.write('gene dictionary for %s generated \n' % (bed_file))
    return gene_dictionary

def gene_name_to_coordinate_conversion(gene_file,gene_dictionary):
    """converts gene name to genomic coordinate in the form of a write out to text file

    Input:
    gene_file: absolute filepath to single column text file with gene names, with one gene name per line
    gene_dictionary: output object of gene_dictionary_from_bed_file()

    Output:
    out: absolute filepath to genomic coordinates for gene names queried by user

    """
    with open(gene_file,'r') as infile:
        wd = os.getcwd()
        out = '%s/gene_to_coordinate.txt' % (wd)
        count,query = 0,0
        with open(out, 'w') as outfile:
            for line in infile:
                clean_line = line.lstrip().rstrip()
                parts = clean_line.split()
                if gene_dictionary.has_key(parts[0]):
                    outfile.write('%s \n' % (gene_dictionary[parts[0]]))
                    count += 1
                    query += 1
                elif re.match(r"^chr\w+:\d+-\d+$", parts[0]) or re.match(r"^chr\d+:\d+-\d+$", parts[0]):
                    outfile.write('%s \n' % (parts[0]))
                    count += 1
                    query += 1
                else:
                    query += 1
                    sys.stdout.write('%s is not an identified element in the dictionary, skipping \n' % parts[0])
                    continue

        sys.stdout.write('gene name to coordinate conversion complete for %s of %s queries \n' % (count,query))
        return out

def within_writeout_csv(guides, output_dir, coordinate, off, annotation_signal):
    """function that write 'within' query output to a CSV file

    Input:
    guides: GuideRNA object which contains sgRNAs within a target coordinate
    out_dir = Output directory to which csv file will be written
    coordinate = Coordinate of region queried. This coordinate will be used to name the output file
    off: True or False value from --off which indicates if all offtarget information should be written to file

    Note:
    The output csv files are written in the excel dialect

    """
    if len(guides) != 0:
        label = 'gRNA_1'
        count = 0
        outfile = ('%s%s%s%s' % (output_dir,'/',coordinate,'_within.csv'))
        if off:
            outfile_offtargets = ('%s%s%s%s' % (output_dir,'/',coordinate,'_offtargets_within.csv'))
            offtarget_csv = open(outfile_offtargets,'wb')
            off_target_writer = csv.writer(offtarget_csv,quoting=csv.QUOTE_MINIMAL, dialect='excel')
            off_target_writer.writerow(['offtarget information for %s gRNAs' % coordinate])
            off_target_writer.writerow(['chromosome','target site start coordinate','target site end coordinate','strand','mismatch',
                                        'intended target site'])

        if annotation_signal == 0:
            final_field = 11
        else:
            final_field = 12

        if type(guides[0]) == list or type(guides[0]) == tuple:
            with open(outfile, 'wb') as output_csv:
                writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
                writer.writerow([coordinate])
                writer.writerow(['chromosome', 'target site start coordinate', 'target site end coordinate', 'gRNA',
                                 'cutting efficiency score', 'cutting specificity score','strand',
                                 'offtargets sum', 'offtargets summary', 'annotation', 'gRNA label'])
                for i in guides:
                    if i[10] == 0:
                        count += 1
                        writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],
                                         i[final_field],label+'__'+str(count)])
                    else:
                        off_target_writer.writerow([i[0],str(i[1]),str(i[2]),i[5],str(i[9]),i[3]])
                if off:
                    offtarget_csv.close()

        else:
            guides = get_bed_lines(guides, offtargets=off)
            with open(outfile, 'wb') as output_csv:
                writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
                writer.writerow([coordinate])
                writer.writerow(['chromosome', 'target site start coordinate', 'target site end coordinate', 'gRNA',
                                 'cutting efficiency score','cutting specificity score','strand',
                                 'offtargets sum', 'offtargets summary', 'annotation', 'gRNA label'])
                for i in guides:
                    if i[10] == 0:
                        count += 1
                        writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],
                                         i[final_field],label+'__'+str(count)])
                    else:
                        off_target_writer.writerow([i[0],str(i[1]),str(i[2]),i[5],str(i[9]),i[3]])
                if off:
                    offtarget_csv.close()

        if off:
            final_outfile = ('%s%s%s%s' % (output_dir, '/', coordinate, '_ontargets_and_offtargets_within.csv'))
            dual_files = [outfile, outfile_offtargets]
            csv_files_as_sheets_in_single_excel(dual_files, final_outfile)
            os.remove(outfile)
            os.remove(outfile_offtargets)

    else:
        sys.stdout.write('no sgRNAs found within %s' % (coordinate))

def flanking_writeout_csv(left_guides, right_guides, output_dir, coordinate, off,annotation_signal):
    """function that write 'flanking' query output to a CSV file

    Input:
    guides_left: GuideRNA object which contains sgRNAs left flanking a target coordinate
    guides_right: GuideRNA object which contains sgRNAs right flanking a target coordiante
    out_dir = Output directory to which csv file will be written
    coordinate = Coordinate of region queried. This coordinate will be used to name the output file
    off: True or False value from --off which indicates if all offtarget information should be written to file

    Note:
    The output csv files are written in the excel dialect

    """
    #convert to bed lines if not already

    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    if len(left_guides) == 0 and len(right_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs are found in either left or right flanking regions: exiting \n')
        sys.exit(1)

    elif len(left_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in left flanking region, writing out only right sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_csv(right_guides, output_dir, ('%s%s' % (coordinate,'_right_flank')),off,annotation_signal)

    elif len(right_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in right flanking region, writing out only left sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_csv(left_guides, output_dir, ('%s%s' % (coordinate,'_left_flank')),off,annotation_signal)

    else:
        if type(left_guides[0]) != list and type(left_guides[0]) != tuple:
            guides_left = get_bed_lines(left_guides, offtargets=False)
        elif type(left_guides[0]) == list or type(left_guides[0]) == tuple:
            guides_left = left_guides
        else:
            sys.stderr.write('ERROR: guides_left object is not GuideRNA class \n')
            sys.exit(1)

        if type(right_guides[0]) != list and type(right_guides[0]) != tuple:
            guides_right = get_bed_lines(right_guides, offtargets=False)
        elif type(right_guides[0]) == list or type(right_guides[0]) == tuple:
            guides_right = right_guides
        else:
            sys.stderr.write('ERROR: guides_right object is not GuideRNA class \n')
            sys.exit(1)

        #make empty line to add to bring short list to same length as long list
        filler = tuple(repeat('NA', len(guides_right[0])))

        #adjust shorter list so that it is same length as longer list
        if len(guides_right) > len(guides_left):
            difference = len(guides_right) - len(guides_left)
            guides_left.extend(repeat(filler, difference))

        elif len(guides_left) > len(guides_right):
            difference = len(guides_left) - len(guides_right)
            guides_right.extend(repeat(filler, difference))

        #write out
        count = 0
        start_coordinate, end_coordinate = int(coordinate.split(':')[1].split('-')[0]), \
                                           int(coordinate.split(':')[1].split('-')[1])
        outfile = output_dir + '/' + coordinate + '_flanking.csv'

        with open(outfile, 'wb') as output_csv:
            writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
            writer.writerow([coordinate])
            writer.writerow(['chromosome', 'target site start coordinate', 'target site end coordinate', 'gRNA', 'cutting efficiency score',
                 'cutting specificity score','strand', 'offtargets sum', 'offtarget summary','left flank distance',
                 'annotation', 'left gRNA label', '', '', 'chromosome', 'target site start coordinate', 'target site end coordinate',
                 'gRNA', 'cutting efficiency score','cutting specificity score','strand', 'offtargets sum', 'offtarget summary',
                 'right flank distance', 'annotation', 'right gRNA label', '',
                 'total flanking distance', 'total deletion size','Vidigal & Ventura s_mU6 oligo'])

            for i, j in zip(guides_left, guides_right):
                label = 'gRNA_1'
                count += 1

                if i[1] and j[1] and i[1] != 'NA' and j[1] != 'NA':
                    left_flank_distance, right_flank_distance = (start_coordinate - int(i[1])), (int(j[1]) - end_coordinate)
                    total_flanking_distance = left_flank_distance + right_flank_distance
                    total_deletion_size = (end_coordinate - start_coordinate) + left_flank_distance + right_flank_distance
                elif i[1] and i[1] != 'NA':
                    left_flank_distance, right_flank_distance = (start_coordinate - int(i[1])), 0
                    total_flanking_distance = left_flank_distance + right_flank_distance
                    total_deletion_size = 'NA'
                    right_flank_distance = 'NA'
                elif j[1] and j[1] != 'NA':
                    left_flank_distance, right_flank_distance = 0, (int(j[1]) - end_coordinate)
                    total_flanking_distance = left_flank_distance + right_flank_distance
                    total_deletion_size = 'NA'
                    left_flank_distance = 'NA'
                else:
                    left_flank_distance, right_flank_distance = None, None
                    total_flanking_distance = None
                    total_deletion_size = None

                #Vidigal & Ventura oligo generation
                if 'GAAGAC' in i[3] or 'GTCTTC' in i[3] or 'GAAGAC' in j[3] or 'GTCTTC' in j[3]:
                    writer.writerow(
                        [i[0], str(i[1]), str(i[2]), i[3].replace('NGG',''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                         str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                         str(j[2]), j[3].replace('NGG',''),str(j[4]),str(j[5]),str(j[6]),j[7],j[8],str(right_flank_distance), j[final_field],
                         label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])
                elif i[3] and j[3] and i[3] != 'NA' and j[3] != 'NA':
                    oligo = ('%s%s%s%s%s' % ('CATGCGAGAAAAGCCTTGTTTGG'.lower(),i[3].replace('NGG',''),'GTTTGGGTCTTCGAGAAGACCTCACCG'.lower()
                                            ,j[3].replace('NGG',''),'GTTTTAGAGCTAGAAATAGC'.lower()))
                    writer.writerow(
                        [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                         str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                         str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                         str(right_flank_distance), j[final_field],
                         label + '_right_' + str(count), '\t', total_flanking_distance,total_deletion_size,oligo])
                else:
                    writer.writerow(
                        [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],
                         str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                         str(j[2]), j[3].replace('NGG', ''), str(j[4]),str(j[5]),str(j[6]), j[7],j[8],
                         str(right_flank_distance), j[final_field],
                         label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])

        if off:

            guides_left = get_bed_lines(left_guides, offtargets=True)
            guides_right = get_bed_lines(right_guides, offtargets=True)

            outfile_offtargets = ('%s%s%s%s' % (output_dir, '/', coordinate, '_offtargets_flanking.csv'))
            offtarget_csv = open(outfile_offtargets, 'wb')
            off_target_writer = csv.writer(offtarget_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
            off_target_writer.writerow(['offtarget information for %s gRNAs' % coordinate])
            off_target_writer.writerow(
                ['chromosome', 'target site start coordinate', 'target site end coordinate', 'strand', 'mismatch',
                 'intended target site'])

            flanking_gRNA = [guides_left,guides_right]
            for z in flanking_gRNA:
                for q in z:
                    if q[10] != 0:
                        off_target_writer.writerow([q[0], str(q[1]), str(q[2]), q[5], str(q[9]), q[3]])

            offtarget_csv.close()

            final_outfile = ('%s%s%s%s' % (output_dir, '/', coordinate, '_ontargets_and_offtargets_flanking.csv'))
            dual_files = [outfile, outfile_offtargets]
            csv_files_as_sheets_in_single_excel(dual_files, final_outfile)
            os.remove(outfile)
            os.remove(outfile_offtargets)

def file_existance(input_file):
    """Ensures a file exists at the stated path

    Input:
    input_file: absolute filepath of the file to be processed

    """
    if os.path.exists(input_file):
        return
    else:
        sys.stderr.write('ERROR: file %s does not exist at state path \n' % (input_file))
        sys.exit(1)

def print_bed(bedlines, target, off, filename=None, header=True):
    """Print guideRNAs and their off-targets to file.

    Args:
    bedlines: list of BED lines such as get_bed_lines() output
    target: input from --target (within) (flanks)
    offtarget tag from args.offtargets --off
    filename: name of file where to print results; if None, print to
              stdout (pipe won't work)
    header: if True, include header
    """
    # util.warn_file_exists(filename)
    if len(bedlines) != 0:
        if type(bedlines[0]) != list and type(bedlines[0]) != tuple:
            bedlines = get_bed_lines(bedlines, offtargets=off)

        if filename:
            f = open(filename, 'w')
        else:
            f = sys.stdout
        if header:
            f.write('# %s\n' % '\t'.join(BED_HEADER))
        for bedline in bedlines:
            f.write('%s\n' % '\t'.join('%s' % x for x in bedline))
        if filename:
            f.close()
    else:
        sys.stdout.write('no sgRNAs found in %s queried region \n' % (target))
        return

def flanking_region_sgRNA_query(database, coordinate, flanking_distance):
    """"function executes query_bam on flanking regions around a coordinate defined by the user

    Input:
    database = filepath to the bam file database produced by processer.py
    coordinate = genomic coordinate inputted by the user. It takes the form chrQ:start-end where Q is an integer of str
    flanking_distance = distance flanking the original coordinate a user desires sgRNAs from, --flankdistance

     Convention:
     query_bam is run with offcoords set to True and onebased set to False by default.

     """
    left_coordinate, right_coordinate, chromosome = int(coordinate.split(':')[1].split('-')[0]), \
                                                    int(coordinate.split(':')[1].split('-')[1]), \
                                                    coordinate.split(':')[0]

    if right_coordinate < left_coordinate:
        sys.stderr.write('ERROR: right coordinate less than left coordinate')
        return

    else:

        if (int(left_coordinate) - int(flanking_distance)) > 0:
            left_flank_coordinate = int(left_coordinate) - int(flanking_distance)
            left_flank_query_region = chromosome + ':' + str(left_flank_coordinate) + '-' + str(left_coordinate)
            left_flank_sgRNAs = query_bam(database, left_flank_query_region, offcoords=True, onebased=False)
            left_flank_sgRNAs = ucsc_compatibility_insurance(left_flank_sgRNAs)
        else:
            sys.stderr.write('WARNING: left flank coordinate not defined, no sgRNAs produced \n')
            left_flank_sgRNAs = []

        right_flank_coordinate = int(right_coordinate) + int(flanking_distance)
        right_flank_query_region = chromosome + ':' + str(right_coordinate) + '-' + str(right_flank_coordinate)
        right_flank_sgRNAs = query_bam(database, right_flank_query_region, offcoords=True, onebased=False)
        right_flank_sgRNAs = ucsc_compatibility_insurance(right_flank_sgRNAs)

        return left_flank_sgRNAs,right_flank_sgRNAs

def csv_files_as_sheets_in_single_excel(dual_files,final_outfile):
    wb = xlwt.Workbook()
    count = 0
    sheet_names = ['ontarget', 'offtarget']
    for filename in dual_files:
        (f_path, f_name) = os.path.split(filename)
        (f_short_name, f_extension) = os.path.splitext(f_name)
        f_short_name,f_extension = sheet_names[count],f_extension
        ws = wb.add_sheet(f_short_name)
        spamReader = csv.reader(open(filename, 'rb'))
        for rowx, row in enumerate(spamReader):
            for colx, value in enumerate(row):
                ws.write(rowx, colx, value)
        count += 1
    wb.save(final_outfile)

def flanking_writeout_csv_bed_input(guides_left,guides_right,outfile,coordinate,off,label,annotation_signal,genomic_coordinate):
    """function that write 'flanking' query output to a CSV file

    Input:
    guides_left: GuideRNA object which contains sgRNAs left flanking a target coordinate
    guides_right: GuideRNA object which contains sgRNAs right flanking a target coordiante
    outfile: absolute filepath to an open file which will accept output.
    coordinate: genomic coordinates of queried target region
    off: True or False value from --off which indicates if all offtarget information should be written to file
    label: label: Unique identifier of each sgRNA pair

    Note:
    The output csv files are written in the excel dialect

    """
    ####################
    #annontation signal#
    ####################
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    #####################################
    #convert to bed lines if not already#
    #####################################
    if len(guides_left) == 0 and len(guides_right) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in flanking regions')
        return

    elif len(guides_left) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in left flanking region, writing out only right sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_csv_bed_input(guides_right,outfile,off,label,annotation_signal,genomic_coordinate)

    elif len(guides_right) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in right flanking region, writing out only left sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_csv_bed_input(guides_left,outfile,off,label,annotation_signal,genomic_coordinate)

    else:
        if type(guides_left[0]) != list and type(guides_left[0]) != tuple:
            guides_left = get_bed_lines(guides_left, offtargets=off)

        if type(guides_right[0]) != list and type(guides_right[0]) != tuple:
            guides_right = get_bed_lines(guides_right, offtargets=off)

        # make empty line to add to bring short list to same length as long list
        filler = tuple(repeat('NA', len(guides_right[0])))

        # adjust shorter list so that it is same length as longer list
        if len(guides_right) > len(guides_left):
            difference = len(guides_right) - len(guides_left)
            guides_left.extend(repeat(filler, difference))

        elif len(guides_left) > len(guides_right):
            difference = len(guides_left) - len(guides_right)
            guides_right.extend(repeat(filler, difference))

        ###########
        #write out#
        ###########
        count = 0
        output_csv = outfile
        start_coordinate,end_coordinate = int(coordinate.split(':')[1].split('-')[0]),\
                                          int(coordinate.split(':')[1].split('-')[1])
        writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
        writer.writerow([genomic_coordinate])
        writer.writerow(['chromosome','target site start coordinate','target site end coordinate','gRNA','cutting efficiency score',
                        'cutting specificity score','strand','offtargets sum','offtarget summary','left flank distance',
                        'annotation','left gRNA label','','','chromosome','target site start coordinate','target site end coordinate',
                        'gRNA','cutting efficiency score','cutting specificity score','strand','offtargets sum','offtarget summary',
                        'right flank distance','annotation','right gRNA label','',
                        'total flanking distance','total deletion size','Vidigal & Ventura s_mU6 oligo'])
        for i, j in zip(guides_left, guides_right):
            count += 1
            if i[1] and j[1] and i[1] != 'NA' and j[1] != 'NA':
                left_flank_distance,right_flank_distance = (start_coordinate-int(i[1])),(int(j[1])-end_coordinate)
                total_flanking_distance = left_flank_distance+right_flank_distance
                total_deletion_size = (end_coordinate-start_coordinate)+left_flank_distance+right_flank_distance
            elif i[1] and i[1] != 'NA':
                left_flank_distance,right_flank_distance=(start_coordinate-int(i[1])),0
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                right_flank_distance = 'NA'
            elif j[1] and j[1] != 'NA':
                left_flank_distance, right_flank_distance = 0, (int(j[1]) - end_coordinate)
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                left_flank_distance = 'NA'
            else:
                left_flank_distance, right_flank_distance = None, None
                total_flanking_distance = None
                total_deletion_size = None

            ####################################
            #Vidigal & Ventura oligo generation#
            ####################################
            if 'GAAGAC' in i[3] or 'GTCTTC' in i[3] or 'GAAGAC' in j[3] or 'GTCTTC' in j[3]:
                writer.writerow(
                    [i[0],str(i[1]),str(i[2]),i[3].replace('NGG', ''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])
            elif i[3] and j[3] and i[3] != 'NA' and j[3] != 'NA':
                oligo = ('%s%s%s%s%s' % (
                'CATGCGAGAAAAGCCTTGTTTGG'.lower(), i[3].replace('NGG', ''), 'GTTTGGGTCTTCGAGAAGACCTCACCG'.lower()
                , j[3].replace('NGG', ''), 'GTTTTAGAGCTAGAAATAGC'.lower()))
                writer.writerow(
                    [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]),j[3].replace('NGG', ''),str(j[4]),str(j[5]),str(j[6]),j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size, oligo])
            else:
                writer.writerow(
                    [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])

        #writer.writerow('\n')

def flanking_writeout_csv_bed_input_offtargets(left_guides,right_guides,outfile,coordinate,off,label,offtarget_csv,annotation_signal,genomic_coordinate,site_guides_left='',site_guides_right=''):
    """function that write 'flanking' query output to a csv file

    Input:
    guides_left = Output of guidequery.query_bam for a left flank region of the target coordinates
    guides_right = Output of guidequery.query_bam for a right flank region of the target coordinates
    out_dir = Output directory to which csv file will be written
    coordinate = Coordinate of region queried. This coordinate will be used to name the output file
    off = offtarget tag from args.offtargets --off

    Note:
    The output csv files are written in the excel dialect
    """
    ##################
    #annotation check#
    ##################
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    #####################################
    #convert to bed lines if not already#
    #####################################
    if len(left_guides) == 0 and len(right_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in flanking regions')
        return

    elif len(left_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in left flanking region, writing out only right sgRNAs "within" '
                         'specified flanking region \n')
        if off:
            within_writeout_csv_bed_input_offtargets(right_guides,outfile,off,label,offtarget_csv,annotation_signal,genomic_coordinate)
        else:
            within_writeout_csv_bed_input(right_guides,outfile,off,label,annotation_signal,genomic_coordinate)

    elif len(right_guides) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in right flanking region, writing out only left sgRNAs "within" '
                         'specified flanking region \n')
        if off:
            within_writeout_csv_bed_input_offtargets(left_guides,outfile,off,label,offtarget_csv,annotation_signal,genomic_coordinate)
        else:
            within_writeout_csv_bed_input(left_guides,outfile,off,label,annotation_signal,genomic_coordinate)

    else:
        if type(left_guides[0]) == list or type(left_guides[0]) == tuple:
            guides_left = left_guides
            print guides_left
        elif type(left_guides[0]) != list and type(left_guides[0]) != tuple:
            guides_left = get_bed_lines(left_guides, offtargets=False)
        else:
            sys.stderr.write('ERROR: guides_left object is not GuideRNA class \n')
            sys.exit(1)

        if type(right_guides[0]) == list or type(right_guides[0]) == tuple:
            guides_right = right_guides
        elif type(right_guides[0]) != list and type(right_guides[0]) != tuple:
            guides_right = get_bed_lines(right_guides, offtargets=False)
        else:
            sys.stderr.write('ERROR: guides_left object is not GuideRNA class \n')
            sys.exit(1)

        ########################################################################
        #make empty line to add to bring short list to same length as long list#
        ########################################################################
        filler = tuple(repeat('NA', len(guides_right[0])))

        ##############################################################
        #adjust shorter list so that it is same length as longer list#
        ##############################################################
        if len(guides_right) > len(guides_left):
            difference = len(guides_right) - len(guides_left)
            guides_left.extend(repeat(filler, difference))

        elif len(guides_left) > len(guides_right):
            difference = len(guides_left) - len(guides_right)
            guides_right.extend(repeat(filler, difference))

        ###########
        #write out#
        ###########
        count = 0
        output_csv = outfile
        start_coordinate,end_coordinate = int(coordinate.split(':')[1].split('-')[0]),\
                                          int(coordinate.split(':')[1].split('-')[1])
        writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
        writer.writerow([genomic_coordinate])
        writer.writerow(['chromosome','target site start coordinate','target site end coordinate','gRNA','cutting efficiency score',
                        'cutting specificity score','strand','offtargets sum','offtarget summary','left flank distance',
                        'annotation','left gRNA label','','','chromosome','target site start coordinate','target site end coordinate',
                        'gRNA','cutting efficiency score','cutting specificity score','strand','offtargets sum','offtarget summary',
                        'right flank distance','annotation','right gRNA label','',
                        'total flanking distance','total deletion size','Vidigal & Ventura s_mU6 oligo'])
        for i, j in zip(guides_left, guides_right):
            count += 1
            if i[1] and j[1] and i[1] != 'NA' and j[1] != 'NA':
                left_flank_distance,right_flank_distance = (start_coordinate-int(i[1])),(int(j[1])-end_coordinate)
                total_flanking_distance = left_flank_distance+right_flank_distance
                total_deletion_size = (end_coordinate-start_coordinate)+left_flank_distance+right_flank_distance
            elif i[1] and i[1] != 'NA':
                left_flank_distance,right_flank_distance=(start_coordinate-int(i[1])),0
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                right_flank_distance = 'NA'
            elif j[1] and j[1] != 'NA':
                left_flank_distance, right_flank_distance = 0, (int(j[1]) - end_coordinate)
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                left_flank_distance = 'NA'
            else:
                left_flank_distance, right_flank_distance = None, None
                total_flanking_distance = None
                total_deletion_size = None

            ####################################
            #Vidigal & Ventura oligo generation#
            ####################################
            if 'GAAGAC' in i[3] or 'GTCTTC' in i[3] or 'GAAGAC' in j[3] or 'GTCTTC' in j[3]:
                writer.writerow(
                    [i[0], str(i[1]),str(i[2]),i[3].replace('NGG', ''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])
            elif i[3] and j[3] and i[3] != 'NA' and j[3] != 'NA':
                oligo = ('%s%s%s%s%s' % (
                'CATGCGAGAAAAGCCTTGTTTGG'.lower(), i[3].replace('NGG', ''), 'GTTTGGGTCTTCGAGAAGACCTCACCG'.lower()
                , j[3].replace('NGG', ''), 'GTTTTAGAGCTAGAAATAGC'.lower()))
                writer.writerow(
                    [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size, oligo])
            else:
                writer.writerow(
                    [i[0], str(i[1]), str(i[2]), i[3].replace('NGG', ''), str(i[4]), str(i[5]), str(i[6]), i[7],i[8],
                     str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                     str(j[2]), j[3].replace('NGG', ''), str(j[4]), str(j[5]), str(j[6]), j[7],j[8],
                     str(right_flank_distance), j[final_field],
                     label + '_right_' + str(count), '\t', total_flanking_distance, total_deletion_size])

        #writer.writerow('\n')

        ################################
        #offtarget detailed information#
        ################################
        if off:
            guides_left = get_bed_lines(site_guides_left, offtargets=True)
            guides_right = get_bed_lines(site_guides_right, offtargets=True)

            off_target_writer = csv.writer(offtarget_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
            off_target_writer.writerow(
                ['chromosome', 'target site start coordinate', 'target site end coordinate', 'strand', 'mismatch',
                 'intended target site'])

            flanking_gRNA = [guides_left,guides_right]
            for z in flanking_gRNA:
                for q in z:
                    if q[10] != 0:
                        off_target_writer.writerow([q[0], str(q[1]), str(q[2]), q[5], str(q[9]), q[3]])

def within_writeout_csv_bed_input(guides,outfile,off,label,annotation_signal,genomic_coordinate):
    """function that write 'within' query output to a CSV file

    Input:
    guides: GuideRNA object which contains sgRNAs within a target coordinate
    outfile: absolute filepath to an open file which will accept output.
    off: True or False value from --off which indicates if all offtarget information should be written to file
    label: Unique identifier of each sgRNA

    Note:
    The output csv files are written in the excel dialect

    """
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    if len(guides) == 0:
        sys.stdout.write('no sgRNAs found within')
        return

    else:
        count = 0
        output_csv = outfile
        writer = csv.writer(output_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
        writer.writerow([genomic_coordinate])
        writer.writerow(['chromosome','target site start coordinate','target site end coordinate','gRNA','cutting efficiency score',
                         'cutting specificity score','strand','offtargets sum','offtargets summary','annotation','gRNA label'])
        if type(guides[0]) == list or type(guides[0]) == tuple:
            for i in guides:
                count += 1
                writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],i[final_field],label+'__'+str(count)])
            #writer.writerow('\n')

        else:
            guides = get_bed_lines(guides, offtargets=off)
            for i in guides:
                count += 1
                writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],i[final_field],label+'__'+str(count)])
            #writer.writerow('\n')

def within_writeout_csv_bed_input_offtargets(guides,outfile,off,label,offtarget_csv,annotation_signal,genomic_coordinate):
    """function that write 'within' query output to a csv file

    Input:
    guides = Output of guidequery.query_bam
    out_dir = Output directory to which csv file will be written
    coordinate = Coordinate of region queried. This coordinate will be used to name the output file
    off = offtarget tag from args.offtargets --off

    Note:
    The output csv files are written in the excel dialect

    """
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    if len(guides) == 0:
        sys.stdout.write('no sgRNAs found within')
        return

    count = 0
    #output_csv = outfile
    off_target_writer = csv.writer(offtarget_csv, quoting=csv.QUOTE_MINIMAL, dialect='excel')
    off_target_writer.writerow(['chromosome', 'start coordinate', 'end coordinate', 'strand', 'mismatch',
                                'intended target site'])

    if type(guides[0]) == list or type(guides[0]) == tuple:
        writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL, dialect='excel')
        writer.writerow([genomic_coordinate])
        writer.writerow(['chromosome', 'target site start coordinate', 'target site end coordinate', 'gRNA',
                         'cutting efficiency score','cutting specificity score','strand',
                         'offtargets sum', 'offtargets summary', 'annotation', 'gRNA label'])
        for i in guides:
            if i[10] == 0:
                count += 1
                writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],i[final_field],label+'__'+str(count)])
            else:
                off_target_writer.writerow([i[0], str(i[1]), str(i[2]), i[5], str(i[10]), i[3]])

        #writer.writerow('\n')

    else:
        guides = get_bed_lines(guides, offtargets=off)
        writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL, dialect='excel')
        writer.writerow([genomic_coordinate])
        writer.writerow(['chromosome', 'target site start coordinate', 'target site end coordinate', 'gRNA',
                         'cutting efficiency score','cutting specificity score','strand',
                         'offtargets sum', 'offtargets summary', 'annotation', 'gRNA label'])
        for i in guides:
            if i[10] == 0:
                count += 1
                writer.writerow([i[0],str(i[1]),str(i[2]),i[3].replace('NGG',''),str(i[4]),str(i[5]),str(i[6]),i[7],i[8],i[final_field],label+'__'+str(count)])
            else:
                off_target_writer.writerow([i[0], str(i[1]), str(i[2]), i[5], str(i[9]), i[3]])

        #writer.writerow('\n')

def flanking_writeout_bed_bed_input(guides_left,guides_right,target,off,outfile,coordinate,label,annotation_signal,header=True):
    """Print guideRNAs and their off-targets to BED file.

    Input:
    guides_left: GuideRNA object which contains sgRNAs left flanking a target coordinate
    guides_right: GuideRNA object which contains sgRNAs right flanking a target coordiante
    target: value from --target, which is a string equaling either "within" or "flanks". In this case target == "flanks"
    off: True or False value from --off which indicates if all offtarget information should be written to file
    outfile: absolute filepath to an open file which will accept output.
    coordinate: genomic coordinates of queried target region
    label: Unique identifier of each sgRNA pair
    header: if True, include header

    """
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    if len(guides_left) == 0 and len(guides_right) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in flanking regions')
        return

    elif len(guides_left) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in left flanking region, writing out only right sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_bed_bed_input(guides_right,target,off,outfile,label,header=True)

    elif len(guides_right) == 0:
        sys.stdout.write('WARNING: no sgRNAs found in right flanking region, writing out only left sgRNAs "within" '
                         'specified flanking region \n')
        within_writeout_bed_bed_input(guides_left,target,off,outfile,label,header=True)

    else:
        if type(guides_left[0]) != list and type(guides_left[0]) != tuple:
            guides_left = get_bed_lines(guides_left, offtargets=off)

        if type(guides_right[0]) != list and type(guides_right[0]) != tuple:
            guides_right = get_bed_lines(guides_right, offtargets=off)

        # make empty line to add to bring short list to same length as long list
        filler = tuple(repeat('', len(guides_right[0])))

        # adjust shorter list so that it is same length as longer list
        if len(guides_right) > len(guides_left):
            difference = len(guides_right) - len(guides_left)
            guides_left.extend(repeat(filler, difference))

        elif len(guides_left) > len(guides_right):
            difference = len(guides_left) - len(guides_right)
            guides_right.extend(repeat(filler, difference))

        # write out
        count = 0
        start_coordinate, end_coordinate = int(coordinate.split(':')[1].split('-')[0]), \
                                           int(coordinate.split(':')[1].split('-')[1])

        if header:
            outfile.write('# %s\n' % '\t'.join(BED_HEADER_BED_INPUT_FLANKS))

        for i, j in zip(guides_left, guides_right):
            count += 1
            if i[1] and j[1]:
                left_flank_distance, right_flank_distance = (start_coordinate - int(i[1])), (int(j[1]) - end_coordinate)
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = (end_coordinate - start_coordinate) + left_flank_distance + right_flank_distance
            elif i[1]:
                left_flank_distance, right_flank_distance = (start_coordinate - int(i[1])), 0
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                right_flank_distance = ''
            elif j[1]:
                left_flank_distance, right_flank_distance = 0, (int(j[1]) - end_coordinate)
                total_flanking_distance = left_flank_distance + right_flank_distance
                total_deletion_size = 'NA'
                left_flank_distance = ''
            else:
                left_flank_distance, right_flank_distance = None, None
                total_flanking_distance = None
                total_deletion_size = None
            outfile.write('%s\t'*34  % (i[0], str(i[1]), str(i[2]), i[3], str(i[4]), str(i[5]), str(i[6]), i[7], str(i[8]), str(i[9]),str(i[10]),
                          str(left_flank_distance), i[final_field], label + '_left_' + str(count), '\t', '\t', j[0], str(j[1]),
                          str(j[2]), j[3],str(j[4]), str(j[5]), str(j[6]), j[7], str(j[8]), str(j[9]),str(i[10]),
                          str(right_flank_distance), j[final_field],label + '_right_' + str(count), '\t',total_flanking_distance,
                          total_deletion_size,'\n'))
        #outfile.write('\n')

def within_writeout_bed_bed_input(guides,target,off,outfile,label,annotation_signal,header=True):
    """Print guideRNAs and their off-targets to BED file.

    Input:
    guides: GuideRNA object which contains sgRNAs within a target coordinate
    target: value from --target, which is a string equaling either "within" or "flanks". In this case target == "within"
    off: True or False value from --off which indicates if all offtarget information should be written to file
    outfile: absolute filepath to an open file which will accept output.
    label: Unique identifier of each sgRNA
    header: if True, include header

    """
    if annotation_signal == 0:
        final_field = 11
    else:
        final_field = 12

    if len(guides) != 0:

        count = 0
        if type(guides[0]) != list and type(guides[0]) != tuple:
            guides = get_bed_lines(guides, offtargets=off)

        if header:
            outfile.write('# %s\n' % '\t'.join(BED_HEADER_BED_INPUT_WITHIN))
        for i in guides:
            count += 1
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                          (i[0],str(i[1]),str(i[2]),i[3],str(i[4]),str(i[5]),str(i[6]),i[7],str(i[8]),str(i[9]),i[final_field],
                           label+'__'+str(count)))

    else:
        sys.stdout.write('no sgRNAs found in %s queried region \n' % (target))
        return

def annotate_bed(bedlines, genome):
    """Annotate BED lines with genomic features.

    Args:
    bedlines: list of BED lines such as get_bed_lines() output
    genome: interval tree instance as returned by util.create_annot_inttree()

    Return:
    same list bedlines as input, but to each tuple a str with annotations added
    """
    if bedlines:
        if type(bedlines[0]) != list and type(bedlines[0]) != tuple:
            bedlines = get_bed_lines(bedlines, offtargets=False)

        for i,t in enumerate(bedlines):
            chrom, start, end = t[0], t[1], t[2]
            annots = None
            inttree = genome.get(chrom)
            if inttree:
                annots = inttree.find(start, end)
            annots = '*' if not annots else ','.join(annots)
            t = t + (annots,)
            bedlines[i] = t
        return bedlines
    else:
        return bedlines

def query_region_check(c):
    """function checks the query input range to ensure the range is not over 1000000 bases.

    Input:
    c = coordinates entered into guidequery through the -c (coords) parameter.

    Convention:
    coordinates take the form chr#:start-end

    Note:
    coordinate range is restricted to 1000000 bases for computational efficiency on the user's machine. Having
    declared this, there is no intrinsic prohibition against query regions beyond 1000000 bases aside from
    user machine considerations.

    """
    if ':' not in c or '-' not in c:
        sys.stderr.write('ERROR: input coordinates %s are not in proper chr#:start-end format \n' % (c))
        return 1
        #sys.exit(1)
    else:
        start,end = int(c.split(':')[1].split('-')[0]),int(c.split(':')[1].split('-')[1])
        if end <= start:
            sys.stderr.write('ERROR: end coordinate %s is less than or equal to start coordinate %s \n' % (end,start))
            return 1
            #sys.exit(1)
        elif (end-start) > 1000000:
            sys.stderr.write('ERROR: query region %s is beyond 1 megabase, please limit query to <= 1000000 bases \n' % (c))
            return 1
            #sys.exit(1)
        else:
            return 0

def bed_labels_unique(bed_file):
    """ensure that labels in BED file are not redundant

    Input:
    bed_file: absolute filepath to a BED file

    Convention:
    The task is achieved by counting the amount of lines in the file and counting the amount of unique 4th field labels.
    if these numbers do not equate then the function returns 1, if they equate the function returns 0. This function
    is only called if input_file_format_identifier returns the signal for BED file. This function is not called for txt
    or GTF files

    """
    cmd1 = ('%s%s' % ('wc -l ', bed_file))
    cmd2 = ('%s%s%s' % ("awk '{print $4}' ", bed_file, " | sort | uniq | wc -l "))

    line_in_file = subprocess.check_output(cmd1, shell=True)
    unique_labels_in_file = subprocess.check_output(cmd2, shell=True)

    lines = int(line_in_file.strip().split()[0])
    labels = int(unique_labels_in_file.strip().split()[0])

    if lines == labels:
        sys.stdout.write('labels in %s are unique \n' % (bed_file))
        return 0
    else:
        sys.stdout.write('labels in %s are not unique, proceed to arbitrary labeling \n' % (bed_file))
        return 1

def input_file_format_identfier(input_file):
    """Verifies input file is an acceptable format that can be processed by batch_query

    Input:
    input_file: absolute filepath of the .bed, .gtf, or .txt file containing only genomic coordinates to be processed

    Convention:
    The function will check if the file is first a txt file. It will look to ensure only only field (column) is present
    with one genomic coordinate per line. The genomic coordinates are expected to take the form chrQ:start-end where
    Q can be either a number of letter. If only one field exists but the data is not in genomic coordinate format then
    the function rejects the file. Otherwise if only one field exists and it has genomic coordinates it is processed.
    If more than one field exists then the first field is checked to ensure it has chr at the start of every string.
    BED files are required to have chromosome information in the first field. GTF files do not require this, but often
    have the first field have chromosome information. If the first field does not have chromosome information then the
    file is rejected. However, if the first field has chromosome information then the file is checked for start and end
    coordinates. In BED files the second and third field respectively have this information. Provided end coordinates
    are always > start coordinates the file will be accepted and processed. In GTF files the fourth and fifth field have
    start and end coordinates respectively. Provided end coordinates are always > start coordinates the file will be
    accepted and processed. If end coordinates !> start coordinates then the file will be rejected. Files do NOT need to
    have .txt, .bed, .gtf extensions to be processed, but every line in the file must conform to .txt, .bed, .gtf format.
    Therefore it is appropiate to strip any headers from the input file before giving it to the package

    """
    with open(input_file) as batch_file:
        batch_file_format = 'reject'
        for line in batch_file:
            if line != '\n':
                clean_line = line.lstrip().rstrip()
                parts = clean_line.split()
                if len(parts) == 1:
                    if re.match(r"^chr\w+:\d+-\d+$", parts[0]) or re.match(r"^chr\d+:\d+-\d+$", parts[0]):
                        batch_file_format = 'txt'
                    else:
                        sys.stderr.write('ERROR: %s is not a permitted file format, ensure only column in .txt file has genomic '
                                         'coordinates or format chr#:start-end \n')
                        batch_file_format = 'reject'
                        break
                        #sys.exit(1)
                else:
                    if 'chr' or 'CHR' in parts[0] and len(parts) >= 3:
                        try:
                            parts[1], parts[2] = int(parts[1]), int(parts[2])
                            if int(parts[2]) - int(parts[1]) > 0:
                                batch_file_format = 'bed'
                            else:
                                sys.stderr.write(
                                    'ERROR: %s is not a permitted file format, ensure 2nd and 3rd fields are start and'
                                    ' end coordinates respectively for a .bed file \n'
                                    % (batch_file))
                                batch_file_format = 'reject'
                                break
                                #sys.exit(1)

                        except ValueError:
                            try:
                                parts[3], parts[4] = int(parts[3]), int(parts[4])
                                if int(parts[4]) - int(parts[3]):
                                    batch_file_format = 'gtf'
                            except IndexError:
                                sys.stderr.write(
                                    'ERROR: %s is not a permitted file format, ensure 3rd and 4th fields are start and'
                                    ' end coordinates respectively for a .gtf file \n' % (input_file))
                                batch_file_format = 'reject'
                                break
                                #sys.exit(1)
                            except:
                                sys.stderr.write('ERROR: %s is not a permitted file format \n' % (input_file))
                                batch_file_format = 'reject'
                                break

                    else:
                        sys.stderr.write('ERROR: %s is not a permitted file format, please enter .bed, .gtf (with first field a '
                                         'chromosome name like chr:#), or .txt file with genomic coordinates' % (input_file))
                        batch_file_format = 'reject'
                        break
                        #sys.exit(1)
            else:
                sys.stderr.write('WARNING: newline character encountered in %s \n' % (input_file))
                continue

    return batch_file_format

def batch_query(input_file, output_dir, bam_file, target, flankdistance, select, sort, n, formating, off,
                annot='', line_limit=100000000,
                bedlines=False, off_file=False):


    """Executes guidequery functions on a .bed, .gtf, or .txt file containing only genomic coordinates

    Input:
    input_file: absolute filepath of the .bed, .gtf, or .txt file containing only genomic coordinates to be processed
    output_dir: output directory for the output file written from the function
    bam_file: GuideScan database BAM file, produced from running processer.py with user options
    target: string having value "within" or "flanks". Determines if sgRNAs are queried for within a region or flanking
    a region
    flankdistance: the distance upstream and downstream of a target coordinate that the software should look for sgRNAs.
            By default this value is 1000
    select: string having value "score" or "offtarget". This option asks the software to choose best n sgRNAs using a
    double sort of either score then offtargets (score) or offtargets then score (offtargets).
    sort: string having value "score", "offtargets", or "coordinates". All sgRNAs from queried regions are sorted by
    one of these metrics
    n: the amount of best sgRNAs that should be selected when the select option is engaged
    formating: the output file of the function can either be a .csv (csv) or a .bed (bed) file. Default is csv.
    off: if true full off-target information is written out to output file, if false it is not.
    annot: absolute filepath to bed file with annotations
    line_limit: determines how many lines the function will process before terminating
    bedlines: if true then bedlines object from get_bed_lines are returned
    off_file: if true then a composite file with detailed offtarget information for sgRNAs is produced: off must be True.

    """

    #############
    # output file#
    #############
    outfile_path = ('%s%s%s' % (output_dir, '/GuideScan_batch_output.', formating))
    outfile = open(outfile_path, 'w')

    ###########################
    # successful region queries#
    ###########################
    within_guides_regions, flanking_guides_region, no_guides_region = 0, 0, 0

    ##################################
    # ensure file format is acceptable#
    ##################################
    file_format = input_file_format_identfier(input_file)
    annotate_dict = defaultdict(list)

    ###############
    # gRNA labeling#
    ###############
    if file_format == 'bed':
        label_signal = bed_labels_unique(input_file)
    elif file_format == 'reject':
        sys.stderr.write('ERROR: %s is not an acceptable file format. Provide .bed, .gtf, or single column .txt file with '
                         'genomic coordinates \n' % (input_file))
        return
    else:
        label_signal = 1

    ######################################################
    # organism chromosome and chromosome length dictionary#
    ######################################################
    verification_dictionary = organism_chromosome_identity_and_length_confirmation(bam_file)

    ##########################
    # annotation interval tree#
    ##########################
    if annot:
        try:
            if annot.keys(): #web-interface modification to allow write-out of annotations
                genome = annot
        except:
            genome = util.create_annot_inttree(annot)
            if genome:
                pass
            else:
                sys.stderr.write('WARNING: file %s not recognized, proceed '
                                 'without annotations \n' % annot)
    #################
    # bedlines return#
    #################
    if target == "within":
        bedlines_within_lst = []
        genomic_coordinate_lst = []
        bedlines_within_OT_lst = []
    elif target == "flanks" or target == "flanking":
        bedlines_left_lst, bedlines_right_lst = [], []
        genomic_coordinate_lst = []
        bedlines_left_OT_lst, bedlines_right_OT_lst = [], []
    else:
        sys.stderr.write('ERROR: %s is not a recognized target parameter: enter "within" or "flanks" \n' % target)
        sys.exit(1)

    if off:
        outfile_offtargets = ('%s%s' % (output_dir, '/GuideScan_offtargets_within.csv'))
        offtarget_csv = open(outfile_offtargets, 'wb')

    ####################
    # processing of data#
    ####################
    with open(input_file, 'r') as infile:
        guide_counter = 1
        line_counter = 0
        for line in infile:  # while loop here to truncate file read-in or if statement
            if line != '\n':
                line_counter += 1
                if line_counter <= line_limit:
                    clean_line = line.lstrip().rstrip()
                    parts = clean_line.split()

                    #######################################################
                    # ensure genomic coordinate belongs to queried organism#
                    #######################################################
                    if file_format == 'bed':
                        if verification_dictionary.has_key(parts[0]):
                            if parts[1] != parts[2]:
                                if int(parts[1]) < int(parts[2]):
                                    if int(parts[1]) <= verification_dictionary[parts[0]] and int(parts[2]) <= \
                                            verification_dictionary[parts[0]]:
                                        genomic_coordinate = '%s:%s-%s' % (parts[0], parts[1], parts[2])
                                    else:
                                        sys.stderr.write('WARNING: %s and/or %s fall off %s which has length %s \n' % (
                                            parts[1], parts[2], parts[0], verification_dictionary[parts[0]]))
                                        continue
                                else:
                                    sys.stderr.write(
                                        'WARNING: start coordinate %s is less than end coordinate %s in %s:%s-%s \n' %
                                        (parts[1], parts[2], parts[0], parts[1], parts[2]))
                                    continue
                            else:
                                sys.stderr.write('WARNING: start coordinate %s equals end coordinate %s in %s:%s-%s \n' %
                                                 (parts[1], parts[2], parts[0], parts[1], parts[2]))
                                continue
                        else:
                            sys.stderr.write('WARNING: %s is not a recognized chromosome for this organism \n' % (parts[0]))
                            continue
                            # genomic_coordinate = '%s:%s-%s' % (parts[0], parts[1], parts[2])

                    elif file_format == 'gtf':
                        if verification_dictionary.has_key(parts[0]):
                            if parts[3] != parts[4]:
                                if int(parts[3]) < int(parts[4]):
                                    if int(parts[3]) <= verification_dictionary[parts[0]] and int(parts[4]) <= \
                                            verification_dictionary[parts[0]]:
                                        genomic_coordinate = '%s:%s-%s' % (parts[0], parts[3], parts[4])
                                    else:
                                        sys.stderr.write('WARNING: %s and/or %s fall off %s which has length %s \n' % (
                                            parts[3], parts[4], parts[0], verification_dictionary[parts[0]]))
                                        continue
                                else:
                                    sys.stderr.write(
                                        'WARNING: start coordinate %s is greater than end coordinate %s in %s:%s-%s \n' %
                                        (parts[3], parts[4], parts[0], parts[3], parts[4]))
                                    continue
                            else:
                                sys.stderr.write('WARNING: start coordinate %s equals end coordinate %s in %s:%s-%s \n' %
                                                 (parts[3], parts[4], parts[0], parts[3], parts[4]))
                                continue
                        else:
                            sys.stderr.write('WARNING: %s is not a recognized chromosome for this organism \n' % (parts[0]))
                            continue
                            # genomic_coordinate = '%s:%s-%s' % (parts[0], parts[3], parts[4])

                    else:  # txt file
                        txt_chrom, txt_start, txt_end = parts[0].split(':')[0], parts[0].split(':')[1].split('-')[0], \
                                                        parts[0].split(':')[1].split('-')[1]
                        if verification_dictionary.has_key(txt_chrom):
                            if txt_start != txt_end:
                                if int(txt_start) < int(txt_end):
                                    if int(txt_start) <= verification_dictionary[txt_chrom] and int(txt_end) <= \
                                            verification_dictionary[txt_chrom]:
                                        genomic_coordinate = parts[0]
                                    else:
                                        sys.stderr.write('WARNING: %s and/or %s fall off %s which has length %s \n' % (
                                            txt_start, txt_end, txt_chrom, verification_dictionary[txt_chrom]))
                                        continue
                                else:
                                    sys.stderr.write(
                                        'WARNING: start coordinate %s is greater than end coordinate %s in %s:%s-%s \n' %
                                        (txt_start, txt_end, txt_chrom, txt_start, txt_end))
                                    continue
                            else:
                                sys.stderr.write('WARNING: start coordinate %s equals end coordinate %s in %s:%s-%s \n' %
                                                 (txt_start, txt_end, txt_chrom, txt_start, txt_end))
                                continue
                        else:
                            sys.stderr.write(
                                'WARNING: %s is not a recognized chromosome for this organism \n' % (txt_chrom))
                            continue
                            # genomic_coordinate = parts[0]

                    ##############################################
                    # ensure query regions are not > 1000000 bases#
                    ##############################################
                    query_region_check_signal = query_region_check(genomic_coordinate)
                    if query_region_check_signal == 0:
                        pass
                    elif query_region_check_signal == 1:
                        continue  # return

                    #####################################################################################
                    # label sgRNA either using 4th field annotation from bed file or arbitrary annotation#
                    #####################################################################################
                    if label_signal == 0:
                        if annotate_dict[genomic_coordinate]:
                            pass
                        else:
                            annotate_dict[genomic_coordinate].append(parts[3])
                    else:
                        if annotate_dict[genomic_coordinate]:
                            pass
                        else:
                            annotate_dict[genomic_coordinate].append('%s_%s' % ('gRNA', guide_counter))
                            guide_counter += 1

                    ##########################
                    # within target sgRNA task#
                    ##########################
                    if target == 'within':
                        guides = query_bam(bam_file, genomic_coordinate, offcoords=off)
                        guides = ucsc_compatibility_insurance(guides)
                        annotation_signal = 0

                        if select:
                            guides = sgRNA_automatic_selection(guides, select, n)
                        else:
                            if sort:
                                guides = sort_guides(guides, sort)

                        if bedlines:
                            bedlines_within = get_bed_lines(guides)
                            bedlines_within_lst.append(bedlines_within)
                            genomic_coordinate_lst.append(genomic_coordinate)
                            if off:
                                bedlines_within_OT = get_bed_lines(guides, offtargets=True)
                                bedlines_within_OT_lst.append(bedlines_within_OT)

                        if annot:
                            if genome:
                                guides = annotate_bed(guides, genome)
                                annotation_signal = 1

                        try:
                            len(guides)
                            if len(guides) == 0:
                                no_guides_region += 1
                            else:
                                within_guides_regions += 1
                        except TypeError:
                            no_guides_region += 1

                        if formating == 'csv':
                            if off:

                                within_writeout_csv_bed_input_offtargets(guides, outfile, off,
                                                                         annotate_dict[genomic_coordinate][0],
                                                                         offtarget_csv, annotation_signal,
                                                                         genomic_coordinate)
                            else:
                                within_writeout_csv_bed_input(guides, outfile, off, annotate_dict[genomic_coordinate][0],
                                                              annotation_signal, genomic_coordinate)

                        elif formating == 'bed':
                            within_writeout_bed_bed_input(guides, target, off, outfile,
                                                          annotate_dict[genomic_coordinate][0], annotation_signal,
                                                          header=True)

                        else:
                            sys.stdout.write('WARNING: %s is not a recogonized format, writing to csv \n' % (formating))
                            if off:
                                within_writeout_csv_bed_input_offtargets(guides, outfile, off,
                                                                         annotate_dict[genomic_coordinate][0],
                                                                         offtarget_csv, annotation_signal,
                                                                         genomic_coordinate)
                            else:
                                within_writeout_csv_bed_input(guides, outfile, off, annotate_dict[genomic_coordinate][0],
                                                              annotation_signal, genomic_coordinate)

                    ############################
                    # flanking target sgRNA task#
                    ############################
                    elif target == 'flanks' or target == 'flanking':
                        left_flank_guides, right_flank_guides = flanking_region_sgRNA_query(bam_file, genomic_coordinate,
                                                                                            flankdistance)
                        annotation_signal = 0

                        if select:
                            guides_left, guides_right = sgRNA_automatic_selection(left_flank_guides, select, n, 'left'), \
                                                        sgRNA_automatic_selection(right_flank_guides, select, n, 'right')
                        else:
                            if sort:
                                guides_left, guides_right = sort_guides(left_flank_guides, sort, 'left'), \
                                                            sort_guides(right_flank_guides, sort, 'right')
                            else:
                                guides_left, guides_right = left_flank_guides, right_flank_guides

                        site_guides_left, site_guides_right = guides_left, guides_right

                        if bedlines:
                            bedlines_left, bedlines_right = get_bed_lines(guides_left), get_bed_lines(guides_right)
                            bedlines_left_lst.append(bedlines_left)
                            bedlines_right_lst.append(bedlines_right)
                            genomic_coordinate_lst.append(genomic_coordinate)
                            if off:
                                bedlines_left_OT, bedlines_right_OT = get_bed_lines(guides_left, offtargets=True), \
                                                                      get_bed_lines(guides_right, offtargets=True)
                                bedlines_left_OT_lst.append(bedlines_left_OT)
                                bedlines_right_OT_lst.append(bedlines_right_OT)

                        if annot:
                            if genome:
                                guides_left, guides_right = annotate_bed(guides_left, genome), annotate_bed(guides_right,
                                                                                                            genome)
                                annotation_signal = 1

                        try:
                            len(guides_left), len(guides_right)
                            if len(guides_left) == 0 or len(guides_right) == 0:
                                no_guides_region += 1
                            else:
                                flanking_guides_region += 1
                        except TypeError:
                            no_guides_region += 1

                        if formating == 'csv':
                            if off:
                                flanking_writeout_csv_bed_input_offtargets(guides_left, guides_right, outfile,
                                                                           genomic_coordinate
                                                                           , off, annotate_dict[genomic_coordinate][0],
                                                                           offtarget_csv, annotation_signal,
                                                                           genomic_coordinate,site_guides_left,site_guides_right)
                            else:
                                flanking_writeout_csv_bed_input(guides_left, guides_right, outfile, genomic_coordinate, off,
                                                                annotate_dict[genomic_coordinate][0], annotation_signal,
                                                                genomic_coordinate)

                        elif formating == 'bed':
                            flanking_writeout_bed_bed_input(guides_left, guides_right, target, off, outfile,
                                                            genomic_coordinate,
                                                            annotate_dict[genomic_coordinate][0], annotation_signal,
                                                            header=True)

                        else:
                            sys.stdout.write('%s is not a is not a recognized format, writing to csv \n' % (formating))
                            if off:
                                flanking_writeout_csv_bed_input_offtargets(guides_left, guides_right, outfile,
                                                                           genomic_coordinate, off,
                                                                           annotate_dict[genomic_coordinate][0],
                                                                           offtarget_csv, annotation_signal,
                                                                           genomic_coordinate)
                            else:
                                flanking_writeout_csv_bed_input(guides_left, guides_right, outfile, genomic_coordinate, off,
                                                                annotate_dict[genomic_coordinate][0], annotation_signal,
                                                                genomic_coordinate)

                else:
                    sys.stdout.write('batch processing limit of %s lines hit \n' % (line_limit))
                    break
            else:
                sys.stderr.write('WARNING: newline character encountered in processing, skipping to next line \n')
                continue

    outfile.close()
    sys.stdout.write('batch processing of %s finished and %s lines processed \n' % (input_file, line_counter))
    if target == 'within':
        sys.stdout.write('%s within queries processed, found gRNAs for %s queried coordinates \n' %
                         (line_counter, within_guides_regions))
    elif target == 'flanks' or target == 'flanking':
        sys.stdout.write(('%s flanking queries processed, found gRNAs in both flanks for %s queried coordinates \n' %
                          (line_counter, flanking_guides_region)))

    ###########################
    # detailed off-target --off#
    ###########################
    if off and off_file:
        offtarget_csv.close()
        final_outfile = ('%s%s' % (output_dir, '/GuideScan_batch_ontargets_and_offtargets.csv'))
        dual_files = [outfile_path, outfile_offtargets]
        try:
            csv_files_as_sheets_in_single_excel(dual_files, final_outfile)
            os.remove(outfile_path)
            os.remove(outfile_offtargets)
        except ValueError:
            sys.stderr.write('WARNING: Unable to merge ontarget and offtarget files; keep files seperate \n')

    #######################
    # website specific code#
    #######################
    if bedlines:
        if target == "within":
            if off:
                sys.stdout.write('bedlines and off options engaged, returning 5 objects')
                return bedlines_within_lst, genomic_coordinate_lst, bedlines_within_OT_lst, line_counter, within_guides_regions
            else:
                sys.stdout.write('bedlines option engaged, returning 4 objects')
                return bedlines_within_lst, genomic_coordinate_lst, line_counter, within_guides_regions
        elif target == "flanks" or target == "flanking":
            if off:
                sys.stdout.write('bedlines and off options engaged, returning 7 objects')
                return bedlines_left_lst, bedlines_right_lst, genomic_coordinate_lst, bedlines_left_OT_lst, \
                       bedlines_right_OT_lst, line_counter, flanking_guides_region
            else:
                sys.stdout.write('bedlines option engaged, returning 5 objects')
                return bedlines_left_lst, bedlines_right_lst, genomic_coordinate_lst, line_counter, flanking_guides_region

def cutting_specificity_verification(bedlines):
    """verifies that cutting specificity scores are present in the queried database

    Input:
    bedlines = output of guidequery.get_bed_lines

    Convention:
    If no specificity scores exist in a database then the set of scores will contain only one element whose identity
    is '*'. This will cause the function to return a 1. If more than one element exists and that element's identity
    is not '*' then 0 will be returned indicating specificity scores exist in the database

    """
    cutting_efficiency_list = []
    for i in bedlines:
        cutting_efficiency_list.append(i[5])

    cutting_efficiency_set = set(cutting_efficiency_list)
    if len(cutting_efficiency_set) == 1 and '*' in cutting_efficiency_set:
        return 1
    else:
        return 0

def cutting_efficiency_verification(bedlines):
    """verifies that Rule 2 Set cutting efficiency scores are present in the queried database

    Input:
    bedlines = output of guidequery.get_bed_lines

    Convention:
    If no Rule 2 Set scores exist in a database then the set of scores will contain only one element whose identity
    is '*'. This will cause the function to return a 1. If more than one element exists and that element's identity
    is not '*' then 0 will be returned indicating Rule 2 Set scores exist in the database

    """
    cutting_efficiency_list = []
    for i in bedlines:
        cutting_efficiency_list.append(i[4])

    cutting_efficiency_set = set(cutting_efficiency_list)
    if len(cutting_efficiency_set) == 1 and '*' in cutting_efficiency_set:
        return 1
    else:
        return 0

def organism_chromosome_identity_and_length_confirmation(bam_file):
    verification_dict = {}
    bamfile = pysam.AlignmentFile(bam_file,'rb')
    chrom_data = bamfile.header
    for i in range(len(chrom_data['SQ'])):
        verification_dict[chrom_data['SQ'][i]['SN']] = chrom_data['SQ'][i]['LN']
    bamfile.close()
    return verification_dict

def query_bam(guidesbam, coords, offcoords=False, onebased=False):
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

    # create dictionary to verify query coordinates correspond to organism chromosome identity and length
    chromosome_dictionary = organism_chromosome_identity_and_length_confirmation(guidesbam)

    if isinstance(coords, basestring):
        coords = coords.strip().split(':')
        coords = (coords[0],) + tuple(int(x) for x in coords[1].split('-'))
    chrom, start, end = coords

    # confirm coordinate exists in genome
    if chromosome_dictionary.has_key(chrom):
        if start <= chromosome_dictionary[chrom] and end <= chromosome_dictionary[chrom]:
            pass
        else:
            sys.stderr.write('WARNING: %s and/or %s fall off %s which has length %s \n' % (
            start, end, chrom, chromosome_dictionary[chrom]))
            return
    else:
        sys.stderr.write('WARNING: %s is not a recognized chromosome for this organism \n' % (chrom))
        return

    if onebased:
        start -= 1  # switch to 0-based coordinates, as expected by pysam
        end -= 1  # and used internally in our pipeline (in the trie etc.)
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
        # ADD CONDITIONALS HERE TO DEAL WITH THESE PARAMETERS ABSENCE
        score = tags.get('ds', -1)
        specificity = tags.get('cs', -1)

        try:
            if util.is_number(score):
                score = float(score)
            else:
                score = -1
        except ValueError:
            sys.stderr.write('WARNING: score %s could not be converted to float, skipping region \n' % (score))
            return

        try:
            if util.is_number(specificity):
                specificity = float(specificity)
            else:
                specificity = -1
        except ValueError:
            sys.stderr.write(
                'WARNING: specificity %s could not be converted to float, skipping region \n' % (specificity))

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
                            offdist=offdist, maxoffcount=maxoffcount, specificity=specificity,
                            offtargets=offtargets)
        guides.append(guiderna)
    bam.close()
    return guides

def get_bed_lines(guides, offtargets=False):
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

def sort_guides(guides,order,flank_direction='center'):
    """function sorts the output of guidequery.query_bam according to user preferences

    Input:
    guides: GuideRNA object which contains sgRNAs. Output of guidequery.query_bam
    order: String equaling 'offtargets', 'score', specificity, or 'coordinate'. When order == 'offtargets' the output is sorted
            according to the amount of enumerated offtargets for each sgRNA in the set. The sgRNAs will be listed in
            ascending order, with the sgRNA with the fewest offtargets appearing first. When order == 'score' the
            output is sorted according to the computed Rule Set 2 on-target cutting efficiency score from Doench et al.
            Nature Biotechnology 2016. The sgRNAs will be listed in descending order, with the largest
            score appearing first. When order == 'coordinates' the output is sorted according to the sgRNA's proximity
            to the target coordinates. The sgRNAs closest to the given coordinate boundary will be listed first. Order,
            by default, is sort by coordinates. When order == 'specificity' the output is sorted according to the computed
            cutting specificity score (CFD) from Doench et al. Nature Biotechnology 2016. The sgRNAs will be listed in
            descending order, with the largest score appearing first.

    Convention:
    offtargets parameter in guidequery.get_bed_lines is set to False

    Note:
    Rule Set 2 on-target cutting efficencey score can be computed and appended to the sgRNA database if the
    database of sgRNAs is composed of 20mers using the Cas9 enzyme

    CFD specificity score can be computed and appended to the sgRNA database if the database of sgRNAs is composed of
    20mers using the Cas9 enzyme

    """
    if len(guides) == 0:
        sys.stdout.write('no sgRNAs found in a target region \n')
        return guides

    bedlines = get_bed_lines(guides, offtargets=False)

    # if sort based on score selected, ensure scores are present
    score_verification = cutting_efficiency_verification(bedlines)
    specificity_validation = cutting_specificity_verification(bedlines)
    if score_verification == 1 and order == 'score':
        sys.stderr.write('WARNING: cutting efficiency scores not in bam file, selection done by coordinates \n')
        return guides
    elif specificity_validation == 1 and order == 'specificity':
        sys.stderr.write('WARNING: cutting specificity scores not in bam file, selection done by coordinates \n')
        return guides
    else:
        # capture data needed to reconstruct GuideRNA class
        offdist = guides[0].offdist
        maxoffcount = guides[0].maxoffcount

        # off-target information captured for each sgRNA before class break
        offtarget_dictionary_GuideRNA = defaultdict(list)
        for j in range(len(guides)):
            key, value = str(guides[j]).split()[3], guides[j].offtargets
            offtarget_dictionary_GuideRNA[key].append(value)

        if order == 'offtargets':
            bedlines.sort(key=lambda x:
            int(x[7]))

        elif order == 'score':
            bedlines.sort(key=lambda x:
            float(x[4] if util.is_number(x[4]) else -1),
                          reverse=True)

        elif order == 'specificity':
            bedlines.sort(key=lambda x:
            float(x[5] if util.is_number(x[5]) else -1),
                          reverse=True)

        elif order == 'coordinates' and flank_direction == 'center' or order == 'coordinates' and flank_direction == 'right':
            pass

        elif order == 'coordinates' and flank_direction == 'left':
            bedlines = bedlines[::-1]

        else:
            sys.stderr.write(
                'ERROR: %s is not a recognized sorting parameter, please enter coordinates, offtargets,'
                ' or score into --sort \n' % (order))
            sys.exit(1)

        # reconstruct GuideRNA class
        bedlines_lst = []
        for i in range(len(bedlines)):
            if type(bedlines[i][4]) != int:
                seq, coord, score, specificity = bedlines[i][3], (
                bedlines[i][0], bedlines[i][1], bedlines[i][2], bedlines[i][6]), -1, bedlines[i][5]
            else:
                seq, coord, score, specificity = bedlines[i][3], (
                bedlines[i][0], bedlines[i][1], bedlines[i][2], bedlines[i][6]), \
                                    (bedlines[i][4] / 100.0), bedlines[i][5]

            for q in offtarget_dictionary_GuideRNA[seq]:
                select_guide = GuideRNA(seq=seq, coord=coord, score=score, specificity=specificity,offdist=offdist, maxoffcount=maxoffcount,
                                        offtargets=q)
                bedlines_lst.append(select_guide)

        return bedlines_lst

def sgRNA_automatic_selection(guides, selection, n, flank_direction='center'):
    """function selects top n sgRNAs for a user according to their preference for offtargets over score of sgRNA
    or vice versa if more than n sgRNAs are returned by guidequery.query_bam

    Input:
    guides: GuideRNA object which contains sgRNAs. Output of guidequery.query_bam
    selection = String equaling 'score' or 'offtargets'. When selection == 'score' the output is first sorted by
    the Rule Set 2 on-target cutting efficiency score from Doench et al. Nature Biotechnology 2016. The sgRNAs with
    the highest scores are sorted to the top of the list. The top n sgRNAs from this sorted list are selected and
    sorted according to their enumerated off-targets. The results of this second sort are returned. This
    selection emphasizes cutting efficiency over potential off-targets. When selection == 'offtargets' the output is first
    sorted by enumerated off-targets. The sgRNAs with the lowest off-targets are sorted to the top of the list. The top
    n sgRNAs from this sorted list are selected and sorted according to their Rule Set 2 on-target cutting efficiency
    score. The results of this second sort are returned. Default is to select according to 'offtargets'. If selection
    == 'coordinates' then the function selects the first n sgRNAs closest to the target site and sorts them by their
    offtargets. If an unknown selection parameter [outside set of (score,offtargets,coordinates)] then the function exits

    Convention:
    offtargets parameter in guidequery.get_bed_lines is set to False within the function.
    This functionality only applies to objects of query_bam with a length >n (more than n sgRNAs outputted per query
    region). For query regions with n sgRNAs or less it is left to the user to choose optimal sgRNAs.

    Note:
    Rule Set 2 on-target cutting efficencey score can be computed and appended to the sgRNA database if the
    database of sgRNAs is composed of 20mers using the Cas9 enzyme

    offtargets are 7th field if specificity scores are present, otherwise they are the 6th field

    """
    # ensure n is a positive integer
    if util.is_number(n) == False:
        sys.stderr.write('%s is not a number, n = 3 \n' % (n))
        n = 3
    elif n <= 0 and util.is_number(n) == True:
        sys.stderr.write('%s is <= 0, n = 3' % (n))
        n = 3
    elif type(n) == float and util.is_number(n) == True:
        sys.stderr.write('%s is a float, converting to integer \n' % (n))
        n = int(n)

    n = int(n)
    if len(guides) > n:

        # capture data needed to reconstruct GuideRNA class
        offdist = guides[0].offdist
        maxoffcount = guides[0].maxoffcount

        # off-target information captured for each sgRNA before class break
        offtarget_dictionary_GuideRNA = defaultdict(list)
        for j in range(len(guides)):
            key, value = str(guides[j]).split()[3], guides[j].offtargets
            offtarget_dictionary_GuideRNA[key].append(value)

        # double sort
        bedlines = get_bed_lines(guides, offtargets=False)

        # if sort based on score selected, ensure scores are present
        score_verification = cutting_efficiency_verification(bedlines)
        specificity_validation = cutting_specificity_verification(bedlines)
        if score_verification == 1 and selection == 'score':
            sys.stderr.write('WARNING: cutting efficiency scores not in bam file, selection done by coordinates on right of query boundry \n')
            selection = 'coordinates'
        elif specificity_validation == 1 and selection == 'specificity':
            sys.stderr.write('WARNING: cutting specificity scores not in bam file, selection done by coordinates on right of query boundry \n')
            selection = 'coordinates'
        else:
            pass

        if selection == 'score':
            bedlines.sort(key=lambda x:
            float(x[4] if util.is_number(x[4]) else -1),
                          reverse=True)
            selects = bedlines[0:n]
            selects.sort(key=lambda x:
            int(x[7] if util.is_number(x[7]) else -1))

        elif selection == 'offtargets':
            bedlines.sort(key=lambda x:
            int(x[7] if util.is_number(x[7]) else -1))
            selects = bedlines[0:n]
            if score_verification == 0:
                selects.sort(key=lambda x:
                float(x[4] if util.is_number(x[4]) else -1),
                             reverse=True)

        elif selection == 'specificity':
            bedlines.sort(key=lambda x:
            float(x[5] if util.is_number(x[5]) else -1),
                          reverse=True)
            selects = bedlines[0:n]
            if score_verification == 0:
                selects.sort(key=lambda x:
                float(x[4] if util.is_number(x[4]) else -1),
                             reverse=True)
            else:
                selects.sort(key=lambda x:
                int(x[7] if util.is_number(x[7]) else -1))

        elif selection == 'coordinates' and flank_direction == 'center' or selection == 'coordinates' and flank_direction == 'right':
            selects = bedlines[0:n]
            selects.sort(key=lambda x:
            int(x[7] if util.is_number(x[7]) else -1))

        elif selection == 'coordinates' and flank_direction == 'left':
            selects = bedlines[(len(bedlines) - n):len(bedlines)]
            selects.sort(key=lambda x:
            int(x[7] if util.is_number(x[7]) else -1))

        else:
            sys.stderr.write(
                'ERROR: %s is not a recognized selection parameter, enter score, offtargets, or coordiantes'
                ' into --select \n' % (selection))
            sys.exit(1)

        # reconstruct GuideRNA class
        selects_lst = []
        for i in range(len(selects)):
            if type(selects[i][4]) != int:
                seq, coord, score, specificity = selects[i][3], (
                    selects[i][0], selects[i][1], selects[i][2], selects[i][6]), -1,selects[i][5]
            else:
                seq, coord, score, specificity = selects[i][3], (
                    selects[i][0], selects[i][1], selects[i][2], selects[i][6]), \
                                    (selects[i][4] / 100.0),selects[i][5]

            for q in offtarget_dictionary_GuideRNA[seq]:
                select_guide = GuideRNA(seq=seq, coord=coord, score=score, specificity=specificity,offdist=offdist,
                                        maxoffcount=maxoffcount,offtargets=q)
                selects_lst.append(select_guide)

        return selects_lst

    else:
        sys.stdout.write('WARNING: <=%s sgRNAs are present in query region, returning all \n' % (n))
        return guides

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

    p = argparse.ArgumentParser(description='Find guideRNAs in precomputed'
                                            ' database stored in BAM format.',
                                epilog='Output is in (extended)'
                                       ' BED format with the following fields:'
                                       + BED_FORMAT_DESCRIPTION +
                                       ' 11: (optionally) list of annotations with genomic'
                                       ' features or "*" if no annotations found.')


    batch_v_single = p.add_mutually_exclusive_group(required=True)
    off_v_annot = p.add_mutually_exclusive_group()
    select_v_sort = p.add_mutually_exclusive_group()
    p.add_argument('-b', dest='bam', required=True,
                   help='path to BAM file with precomputed guideRNAs. REQUIRED', )
    batch_v_single.add_argument('-c', dest='coords',
                                help='coordinates in the form "<chr>:<start>-<end>" ;'
                                     'example: chrX:3364088-3372035, mutually exclusive with --batch')
    batch_v_single.add_argument('--batch', dest='batch_mode',
                                help='absolute filepath to a BED file, GTF file where the first field '
                                     ' (column) contains chromosome information (chr#:), or a txt file '
                                     ' composed of a single field (column) of genomic coordinates of the'
                                     ' following format: chr#:start-end, mutually exclusive with -c')
    batch_v_single.add_argument('--sequence',dest='sequence_file',help='sequence file is fasta file format. '
                                                                       'sequences will be processed through locally '
                                                                       'installed blat binary and sequences with perfect '
                                                                       'matches to a specified genome will be processed. '
                                                                       'If a perfect match does not exist, the individual query '
                                                                       'will not be processed ')
    p.add_argument('--target', default='within', required=True, dest='target',
                   help='get sgRNAs within the target coordinates (within) or flanking the target coordinates'
                        '(flanks) by a distance detailed with the flankdistance parameter. Default is within.')
    p.add_argument('--flankdistance', default=1000, dest='flankdistance',
                   help='the distance flanking both the upstream and downstream regions of a target coordinate.'
                        'Default is 1000.')
    p.add_argument('--one', dest='onebased', action='store_true',
                   help='whether input coordinates are 1-based'
                        ' (default is 0-based); output is always in BED'
                        ' format and 0-based')
    p.add_argument('-o', dest='output', required=True,
                   help='name of output directory. REQUIRED')
    p.add_argument('--output_format', dest='format', default='csv',
                   help='file format for output can either be bed (bed) or csv (csv) format. Default is csv.')
    p.add_argument('--header', action='store_true',
                   help='whether header should be included in output')
    off_v_annot.add_argument('--off', dest='offtargets', action='store_true',
                             help='whether detailed info about off-targets of each'
                                  ' guideRNA should be included in output'
                                  ' (each off-target in a separate line); default is'
                                  ' to include only summary about all off-targets'
                                  ' of a guideRNA. Must be utilized with csv output format (which is default)')
    select_v_sort.add_argument('--sort', dest='sort',
                               help='sort sgRNAs by fewest off-targets (offtargets), highest Rule 2 Set cutting efficiency'
                                    'score (score), highest CFD cutting specificity (specificity), or sgRNAs closest to the target region (coordinates). Default is '
                                    'coordinates. Mutually exclusive with select.')
    select_v_sort.add_argument('--select', dest='select',
                               help='guidequery chooses n optimal sgRNAs based on emphasizing fewest off-targets then '
                                    'sorts by Rule 2 Set cutting efficiency (offtargets) or it chooses n optimal sgRNAs'
                                    'based on emphasizing highest Rule 2 Set cutting efficiency score then sorts '
                                    'by fewest off-targets (score). Also can choose n sgRNAs closes to target and sorts '
                                    'by offtargets (coordinates). Also can choose n sgRNAs with highest CFD score and '
                                    'sorts by offtargets. Mutually exclusive with sort.')
    p.add_argument('-n', dest='n', default=3,
                   help='amount of optimal sgRNAs desired from the --select parameter. Ignored if --select not used. '
                        'Default is 3.')
    off_v_annot.add_argument('--annot', dest='annotfile',
                             help='path to BED file with coordinates of genomic features'
                                  ' that should be used for annotation'
                                  ' (format (tabulated): chrom, start, end, name);'
                                  ' for example, use Table Browser'
                                  ' https://genome.ucsc.edu/cgi-bin/hgTables'
                                  ' to create such BED files of various kinds;'
                                  ' alternatively, use here short names for'
                                  ' preinstalled exon annotations: "hg38" for human,'
                                  ' "dm6" for fly, "mm10" for mouse,'
                                  ' "sacSer3" for yeast, ce11 for c. elegans')
    p.add_argument('--feature_bed_file',dest='features',
                   help='path to BED file with genomic features of interest such that user can upload a BED file with '
                        'feature names (in 4th field) and GuideScan will map the genomic coordinates to the feature ')
    p.add_argument('--fasta_file',dest='fasta',help='fasta file with index located in same directory so that --sequence '
                                                    'queries can be processed')
    p.add_argument('--blat',dest='blat',help='absolute filepath to blat tool. Binary packages for locally running blat '
                                             'are found within GuideScan software at ./blat_binaries')

    args = p.parse_args()

    if args.sequence_file:

        file_existance(args.sequence_file)

        if args.fasta:

            fasta_index_check = file_existance_fasta_index(args.fasta)
            if fasta_index_check == 1:
                sys.stderr.write('WARNING: %s file and/or its index not identified directory\n' % args.fasta)
                return

            else:
                pass

            if args.blat:

                coordinate_file = blat_processing(args.blat, args.fasta, args.sequence_file, args.output, verbose=True)

                if coordinate_file:

                    if args.offtargets:
                        batch_query(coordinate_file, args.output, args.bam, args.target, args.flankdistance, args.select,
                                    args.sort,
                                    args.n,
                                    args.format, args.offtargets, args.annotfile, off_file=True)
                    else:
                        batch_query(coordinate_file, args.output, args.bam, args.target, args.flankdistance, args.select,
                                    args.sort,
                                    args.n, args.format, args.offtargets, args.annotfile, off_file=False)

                    os.remove(coordinate_file)

                else:
                    sys.stderr.write('WARNING: blat is not locally installed; cannot process file\n')

            else:
                sys.stderr.write('WARNING: blat binary is not specified\n')
                return
        else:
            sys.stderr.write('WARNING: fasta file is not specified\n')
            return

    elif args.batch_mode:

        file_existance(args.batch_mode)

        if args.features:
            file_existance(args.features)
            features_dictionary = gene_dictionary_from_bed_file(args.features)
            coordinate_file = gene_name_to_coordinate_conversion(args.batch_mode, features_dictionary)

            if args.offtargets:
                batch_query(coordinate_file, args.output, args.bam, args.target, args.flankdistance, args.select,
                            args.sort,
                            args.n,
                            args.format, args.offtargets, args.annotfile, off_file=True)
            else:
                batch_query(coordinate_file, args.output, args.bam, args.target, args.flankdistance, args.select,
                            args.sort,
                            args.n, args.format, args.offtargets, args.annotfile, off_file=False)

            os.remove(coordinate_file)

        else:
            if args.offtargets:
                batch_query(args.batch_mode, args.output, args.bam, args.target, args.flankdistance, args.select, args.sort,
                            args.n,
                            args.format, args.offtargets, args.annotfile, off_file=True)
            else:
                batch_query(args.batch_mode, args.output, args.bam, args.target, args.flankdistance, args.select, args.sort,
                            args.n, args.format, args.offtargets, args.annotfile, off_file=False)
    else:
        annotation_signal = 0
        query_region_check_signal = query_region_check(args.coords)
        if query_region_check_signal == 0:
            pass
        elif query_region_check_signal == 1:
            return

        if args.target == 'within':

            guides = query_bam(args.bam, args.coords, offcoords=args.offtargets, onebased=args.onebased)
            guides = ucsc_compatibility_insurance(guides)

            try:
                len(guides)
            except TypeError:
                sys.stderr.write('query coordinate is not searchable in queried genome \n')
                return

            if args.select:
                bedlines = sgRNA_automatic_selection(guides, args.select, args.n)
            else:
                if args.sort:
                    bedlines = sort_guides(guides, args.sort)
                else:
                    bedlines = get_bed_lines(guides, offtargets=args.offtargets)

            if args.annotfile:
                genome = util.create_annot_inttree(args.annotfile)
                if genome:
                    bedlines = annotate_bed(bedlines, genome)
                    annotation_signal = 1
                else:
                    sys.stderr.write('WARNING: file %s not recognized, proceed '
                                     'without annotations \n' % args.annotfile)

            if args.format == 'bed':
                if args.select:
                    filename = ('%s%s%s%s' % (args.output, '/GuideScan_Selected_', args.coords, '_within.bed'))
                    print_bed(bedlines, 'within', args.offtargets, filename=filename, header=True)
                else:
                    filename = ('%s%s%s%s' % (args.output, '/', args.coords, '_within.bed'))
                    print_bed(bedlines, 'within', args.offtargets, filename=filename, header=True)
            elif args.format == 'csv':
                if args.select:
                    within_writeout_csv(bedlines, args.output, ('%s%s' % ('GuideScan_Selected_', args.coords)),
                                        args.offtargets, annotation_signal)
                else:
                    within_writeout_csv(bedlines, args.output, args.coords, args.offtargets, annotation_signal)
            else:
                sys.stderr.write('format %s not recognized, writing to csv output \n' % (args.format))
                if args.select:
                    within_writeout_csv(bedlines, args.output, ('%s%s' % ('GuideScan_Selected_', args.coords)),
                                        args.offtargets, annotation_signal)
                else:
                    within_writeout_csv(bedlines, args.output, args.coords, args.offtargets, annotation_signal)

        elif args.target == 'flanks' or args.target == 'flanking':

            left_flank_sgRNAs, right_flank_sgRNAs = flanking_region_sgRNA_query(args.bam, args.coords, args.flankdistance)

            try:
                len(left_flank_sgRNAs)
                len(right_flank_sgRNAs)
            except TypeError:
                sys.stderr.write('query coordinate is not searchable in queried genome \n')
                return

            if args.select:
                bedlines_left = sgRNA_automatic_selection(left_flank_sgRNAs, args.select, args.n, 'left')
                bedlines_right = sgRNA_automatic_selection(right_flank_sgRNAs, args.select, args.n, 'right')
            else:
                if args.sort:
                    bedlines_left, bedlines_right = sort_guides(left_flank_sgRNAs, args.sort, 'left'), \
                                                    sort_guides(right_flank_sgRNAs, args.sort, 'right')
                else:
                    bedlines_left, bedlines_right = get_bed_lines(left_flank_sgRNAs, offtargets=args.offtargets), \
                                                    get_bed_lines(right_flank_sgRNAs, offtargets=args.offtargets)

            if args.annotfile:
                genome = util.create_annot_inttree(args.annotfile)
                if genome:
                    bedlines_left, bedlines_right = annotate_bed(bedlines_left, genome), annotate_bed(bedlines_right,
                                                                                                      genome)
                    annotation_signal = 1
                else:
                    sys.stderr.write('WARNING: file %s not recognized, proceed '
                                     'without annotations \n' % args.annotfile)

            if args.format == 'bed':
                if args.select:
                    filename_left = ('%s%s%s%s' % (args.output, '/GuideScan_Selected_', args.coords, '_left_flank.bed'))
                    filename_right = ('%s%s%s%s' % (args.output, '/GuideScan_Selected_', args.coords, '_right_flank.bed'))
                    print_bed(bedlines_left, 'left_flank', args.offtargets, filename=filename_left, header=True)
                    print_bed(bedlines_right, 'right_flank', args.offtargets, filename=filename_right, header=True)
                else:
                    filename_left = ('%s%s%s%s' % (args.output, '/', args.coords, '_left_flank.bed'))
                    filename_right = ('%s%s%s%s' % (args.output, '/', args.coords, '_right_flank.bed'))
                    print_bed(bedlines_left, 'left_flank', args.offtargets, filename=filename_left, header=True)
                    print_bed(bedlines_right, 'right_flank', args.offtargets, filename=filename_right, header=True)
            elif args.format == 'csv':
                if args.select:
                    flanking_writeout_csv(bedlines_left, bedlines_right, args.output,
                                          ('%s%s' % ('GuideScan_Selected_', args.coords)), args.offtargets,
                                          annotation_signal)
                else:
                    flanking_writeout_csv(bedlines_left, bedlines_right, args.output, args.coords, args.offtargets,
                                          annotation_signal)
            else:
                sys.stderr.write('format %s not recognized, writing to csv output \n' % (args.format))
                if args.select:
                    flanking_writeout_csv(bedlines_left, bedlines_right, args.output,
                                          ('%s%s' % ('GuideScan_Selected_', args.coords)), args.offtargets,
                                          annotation_signal)
                else:
                    flanking_writeout_csv(bedlines_left, bedlines_right, args.output, args.coords, args.offtargets,
                                          annotation_signal)

        else:
            sys.stderr.write('WARNING: target %s not recognized, enter "within" or "flanks" for --target \n'
                             % args.target)

if __name__ == '__main__':
    main()