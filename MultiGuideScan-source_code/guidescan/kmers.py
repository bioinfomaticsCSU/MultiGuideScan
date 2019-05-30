__author__ = 'Tao Li & Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Extract k-mers from genomic sequence and do initial processing.

extract_process_kmers(): main function, call other functions and print log
load_fasta(): load FASTA files
extract_kmers(): extract k-mers and their coordinates from loaded FASTA
print_stats_kmers(): auxiliary function, print statistics of extracted k-mers
shuffle_kmers(): shuffle extracted k-mers, important for guideRNA analysis
sort_count_kmers(): count occurrences in genome of all extracted k-mers
"""

#################
#               #
#   Libraries   #
#               #
#################

import os
import sys
import gzip
import tempfile
import itertools
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq

import util
from multiprocessing import Pool
#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def load_fasta(path):
    """Load FASTA files.

    Args:
    path: path to fasta file or folder with fasta files;
          if folder, will use all .fa and .fasta files found in the folder

    Return:
    iterator over Bio.SeqRecord.SeqRecord objects containing chromosomes
    """
    if os.path.isdir(path):
        filenames = os.listdir(path)
        filenames = ['%s/%s' % (path, f) for f in filenames
                     if f.endswith('.fasta') or f.endswith('.fa')
                        or f.endswith('.fasta.gz')
                        or f.endswith('.fa.gz')]
    elif os.path.isfile(path):
        filenames = [path]
    else:
        raise util.iGuideError('FASTA not found at %s' % path)
    files = [gzip.open(f) if f.endswith('.gz') else open(f)
             for f in filenames]
    fasta = [SeqIO.parse(f, 'fasta') for f in files]
    return itertools.chain(*fasta)

def extract_kmers(name, fasta, length, pams, pampos, filename, chroms=[],
                  minchrlen=10000, processes=1):
    """Extract candidate k-mer guideRNAs with their coordinates from FASTA.

    Convention: coordinate reported is for start position of the whole probe
                in the genome in 0-based coordinates; probe includes guideRNA
                and PAM (PAM can be before or after the guide); then for plus
                strand probe continues to the right, for minus strand probe
                continues to the left

    Args:
    name: project name, a folder with this name containing intermediate and 
            final files in it
    fasta: iterator over Bio.SeqRecord.SeqRecord objects containing chromosomes
           as returned by load_fasta()
    length: length of guideRNAs (not including PAM sequence)
    pams: list of primary and alternative PAM sequences
    pampos: position of PAM ('start' or 'end')
    filename: all k-mers will be written in this file in the format
              '<k-mer followed by PAM> <coordinates>''
    chroms: if not empty, inlcude in analysis only chromosomes with names
            from this list
    minchrlen: include in analysis only chromosomes not shorter than this
    processes: how many processes to use; do not specify more than you have 
                on your system

    Return:
    genome info in the format [(<chromosome name>, <chromosome length>)]
        for all processed chromosomes in the order of processing
    """
    if not pampos in ['start', 'end']:
        raise util.iGuideError("'pampos' argument should be 'start' or 'end'")
    pams_extend = [(pam, pam_seq) for pam in pams
                                  for pam_seq in util.expand_dna_n(pam)]
    pams_extend_rev = [(pam, str(Seq(pam_seq).reverse_complement()))
                       for pam in pams
                       for pam_seq in util.expand_dna_n(pam)]
    genome = []
    fasta_temp = []
    for chrom in fasta:
        if len(chrom) < minchrlen:
            continue
        if chroms and chrom.id not in chroms:
            continue
        genome.append((chrom.id, len(chrom)))
        fasta_temp.append(chrom)

    #Parallelize extracting kmers from the reference genome sequences by user defined 
    #processes or the number of reference sequences
    parts = len(fasta_temp)
    if processes > parts:
        processes = parts

    kmersfiles_temp = [tempfile.NamedTemporaryFile(dir=name,
                                               suffix='.temp%s' % i)
                        for i in range(parts)]

    pool = Pool(processes)
    util.print_log('poolSize %s...' % processes)
    
    for i in range(parts):
        pool.apply_async(extract_kmers_pool,(fasta_temp[i], length, pampos, pams_extend,
                                             pams_extend_rev, kmersfiles_temp[i].name))
    util.print_log('Waiting for all subprocesses done...')
    pool.close()
    pool.join()
    util.print_log('all chromosomes processed')
    
    util.print_log('done, merge all kmers...')
    total_count = 0
    util.warn_file_exists(filename)
    f = gzip.open(filename, 'w')
    for i in range(parts):
        for line in kmersfiles_temp[i]:
            f.write(line)
            total_count += 1

    for file in kmersfiles_temp:
        file.close()
    f.close()
    util.print_log('total k-mers written: %s' % total_count)

    return genome
def extract_kmers_pool(chrom, length, pampos, pams_extend, pams_extend_rev, filename):
    
    f = open(filename, 'w')
    
    util.print_log('process %s' % str((chrom.id, len(chrom))))
    count = 0
    chrom_seq = chrom.seq.upper()
    # plus strand
    for pam, pam_seq in pams_extend:
        pos = -1
        while True:
            pos = chrom_seq.find(pam_seq, pos + 1)
            if pos > -1:
                if pampos == 'start':
                    start = pos + len(pam)
                    coord = pos
                    seq = str(chrom_seq[start:start + length]).upper()
                    if len(seq) == length:  # for boundary cases
                        f.write('%s%s\t%s:%s:+\n' % (pam, seq, chrom.id,
                                                     coord))
                        count += 1
                        total_count += 1
                elif pampos == 'end':
                    start = pos - length
                    coord = start
                    seq = str(chrom_seq[start:start + length]).upper()
                    if len(seq) == length:  # for boundary cases
                        f.write('%s%s\t%s:%s:+\n' % (seq, pam, chrom.id,
                                                     coord))
                        count += 1
            else:
                break
    # minus strand
    for pam, pam_seq in pams_extend_rev:
        pos = -1
        while True:
            pos = chrom_seq.find(pam_seq, pos + 1)
            if pos > -1:
                if pampos == 'start':
                    start = pos - length
                    coord = pos + len(pam) - 1
                    seq = chrom_seq[start:start + length].reverse_complement()
                    seq = str(seq).upper()
                    if len(seq) == length:  # for boundary cases
                        f.write('%s%s\t%s:%s:-\n' % (pam, seq, chrom.id,
                                                     coord))
                        count += 1
                elif pampos == 'end':
                    start = pos + len(pam)
                    coord = pos + len(pam) + length - 1
                    seq = chrom_seq[start:start + length].reverse_complement()
                    seq = str(seq).upper()
                    if len(seq) == length:  # for boundary cases
                        f.write('%s%s\t%s:%s:-\n' % (seq, pam, chrom.id,
                                                     coord))
                        count += 1
            else:
                break
    util.print_log('%s k-mers' % count)
    f.close()

def print_stats_kmers(filename, gnupath=''):
    """Print count of k-mers per chrosomosome per strand.

    Relies on utilities 'gzip', 'cut', 'awk'

    Args:
    filename: name of file with k-mers, assume file is gzipped
    gnupath: path to gnu binary for cut
    """
    gnupath = gnupath.strip()
    if gnupath:
        gnupath = gnupath + '/'
    else:
        gnupath = ''
    f = tempfile.NamedTemporaryFile(dir='.')
    extract_command = "gzip -cd %s | %scut -f2 | " \
                      "awk -F ':' '{print $1\":\"$3}' >%s" \
                      % (filename, gnupath, f.name)
    sys.stdout.write(extract_command)
    os.system(extract_command)
    count = Counter(line.strip() for line in open(f.name))
    print sorted(count.iteritems())
    f.close()

def sort_count_kmers(fileinput, fileoutput, mincount=10, gnupath=''):
    """Sort and count k-mers.

    Relies on utilities 'gzip', 'cut', 'sort', 'uniq', 'awk' available
    in the system.

    Args:
    fileinput: name of file with k-mers, assume file is gzipped
    fileoutput: name of output, will be gzipped
    mincount: show in output only k-mers with at least these many occurrences
    gnupath: path to gnu binaries for cut, sort, uniq; if empty, use defaults;
             using most recent gnu coreutils version can help speed up

    Output:
    produce file with k-mers and counts of how many times each k-mer appeared
    """
    gnupath = gnupath.strip()
    if gnupath:
        gnupath = gnupath + '/'
    else:
        gnupath = ''
    f = tempfile.NamedTemporaryFile(dir='.')
    sort_command = 'gzip -cd %s | %scut -f1 | %ssort >%s' % \
                   (fileinput, gnupath, gnupath, f.name)
    sys.stdout.write(sort_command)
    os.system(sort_command)
    count_command = "%suniq -c %s | awk '{if ($1 >= %s) print $2\" \"$1}'" \
                    " | gzip >%s" % \
                    (gnupath, f.name, mincount, fileoutput)
    sys.stdout.write(count_command)
    os.system(count_command)
    f.close()

def shuffle_kmers(fileinput, fileoutput, gnupath=""):
    """Reproducibly shuffle lines in file.

    Relies on utilities 'gzip' and 'shuf' or 'gshuf' available in the system. Reproducible,
    always produces the same result. For reproducibility uses bytes from
    the input file itself as random source (assume it is random enough).

    Args:
    fileinput: name of file, assume file is gzipped
    fileoutput: name of output, will be gzipped
    gnupath: path to gnu binary for shuf; if empty, use defaults;
             (note: shuf is not available by default on some systems)

    Output:
    produce file with all input lines in randomly shuffled order
    """
    gnupath = gnupath.strip()
    if gnupath:
        gnupath = gnupath + '/'
    else:
        gnupath = ''
    shuf_cmd_check = 'shuf --help'
    gshuf_cmd_check = 'gshuf --help'
    shuf_exist,gshuf_exist = os.system(shuf_cmd_check),os.system(gshuf_cmd_check)
    if shuf_exist == 0:
        shuf_command = 'gzip -cd %s | %sshuf --random-source=%s | gzip >%s' \
                       % (fileinput, gnupath, fileinput, fileoutput)
    elif gshuf_exist == 0:
        shuf_command = 'gzip -cd %s | %sgshuf --random-source=%s | gzip >%s' \
                       % (fileinput, gnupath, fileinput, fileoutput)
    else:
        sys.stderr.write('ERROR: niether shuf nor gshuf are installed on OS, please install: brew install shuf \n')
        sys.exit(1)
    sys.stdout.write(shuf_command)
    os.system(shuf_command)

def extract_process_kmers(name):
    """Extract k-mers from genomic sequence and run initial processing.

    Load project arguments and produce three files:
    extract k-mers from the genome: <name>/<name>_kmers.txt.gz
    shuffle all extracted k-mers: <name>/<name>_kmers_shuffled.txt.gz
    count occurrences of k-mers: <name>/<name>_kmers_counts.txt.gz

    Args:
    name: project name, used to get project args and in all output
    """
    util.print_log('start extract_process_kmers()')
    util.print_log('load arguments...')
    args = util.load_args(name)
    util.print_args(args)
    util.print_log('done')

    util.print_log('load FASTA...')
    util.print_log('load from %s' % args['fasta'])
    fasta = load_fasta(args['fasta'])
    util.print_log('done')

    util.print_log('extract k-mers...')
    kmers_filename = '%s/%s_kmers.txt.gz' % (name, name)
    allpams = [args['pam']] + args['altpam']
    util.print_log('write in file %s' % kmers_filename)
    genome = extract_kmers(name=name, fasta=fasta, length=args['length'],
                           pams=allpams, pampos=args['pampos'],
                           filename=kmers_filename, chroms=args['chrom'],
                           minchrlen=args['minchrlen'], processes=args['processes'])
    sys.stdout.write('genome: %s' % genome)
    util.print_log('save genome info')
    args['genome'] = genome
    util.save_args(args)
    util.print_log('calculate k-mer statistics')
    print_stats_kmers(kmers_filename, gnupath=args['gnupath'])
    util.print_log('done')

    util.print_log('shuffle k-mers...')
    kmers_shuffled_filename = '%s/%s_kmers_shuffled.txt.gz' % (name, name)
    util.print_log('write in file %s' % kmers_shuffled_filename)
    shuffle_kmers(fileinput=kmers_filename, fileoutput=kmers_shuffled_filename,
                  gnupath=args['gnupath'])
    util.print_log('done')

    util.print_log('count k-mers...')
    count_filename = '%s/%s_kmers_counts.txt.gz' % (name, name)
    util.print_log('write in file %s' % count_filename)
    sort_count_kmers(fileinput=kmers_filename, fileoutput=count_filename,
                     mincount=args['maxoffpos'], gnupath=args['gnupath'])
    util.print_log('done')
    return True
