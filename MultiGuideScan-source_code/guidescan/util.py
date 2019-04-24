"""Utilities."""

import os
import sys
import operator
import cPickle
import itertools
import binascii
from datetime import datetime
from pkg_resources import resource_exists, resource_filename

import psutil
import numpy as np
from bx.intervals.intersection import IntervalTree


class iGuideError(Exception):
    pass

class iGuideSAMError(iGuideError):
    pass

def check_file_exists(filename):
    if not os.path.isfile(filename):
        raise iGuideError('file %s is not found' % filename)

def warn_file_exists(filename):
    if os.path.isfile(filename):
        sys.stderr.write('Warning: overwrite existing file %s' % filename)

def is_number(s):
    """Check if input str can be coerced to number"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def print_log(s=''):
    """Print log str with current time and memory usage.

    Args:
    s: str
    """
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 2 ** 20
    sys.stderr.write('[%s %sMb]: %s \n' % (datetime.now(), mem, s))
    sys.stdout.flush()

def print_log2(s=''):
    """Print log str with current time and memory usage.

    Args:
    s: str
    """
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 2 ** 20
    sys.stdout.write('[%s %sMb]: %s \n' % (datetime.now(), mem, s))
    sys.stdout.flush()

def dump(object, filename):
    """Use cPickle.dump to save the object into file 'filename'."""
    f = open(filename, 'wb')
    cPickle.dump(object, f, cPickle.HIGHEST_PROTOCOL)
    f.close()

def load(filename):
    """Use cPickle.load to load the object into file 'filename'."""
    f = open(filename, 'rb')
    object = cPickle.load(f)
    f.close()
    return object

def save_args(args_dict, name=None, warning=True):
    """Save pipeline arguments in file.

    Args:
    args_dict: dict with arguments
    name: to use as dirname and in filename;
          if None or empty, use args_dict['name']
    warning: if True, warn when overwriting existing file

    Pickle in file <name>/<name>_args.pic
    """
    name = name or args_dict['name']
    if not os.path.isdir(name):
        os.makedirs(name)
    filename = '%s/%s_args.pic' % ((name,) * 2)
    if warning:
        warn_file_exists(filename)
    dump(args_dict, filename)

def load_args(name):
    """Load pipeline arguments and return as dict."""
    filename = '%s/%s_args.pic' % ((name,) * 2)
    if not os.path.isfile(filename):
        raise iGuideError('file %s is not found \n' % filename)
    args_dict = load(filename)
    return args_dict

def print_args(args):
    """Print dict of arguments."""
    for k,v in args.iteritems():
        sys.stderr.write('%s : %s \n' % (k, str(v)))

def hamming(s1, s2):
    """Return Hamming distance between two strings of equal length."""
    return sum(map(operator.ne, s1, s2))

def expand_dna_n(s):
    """Replace each occurrence of N in a string by each of ACGT.

    Args:
    s: string composed of letters A, C, G, T, N

    Return a set of strings with all combinations of substitutions
    """
    letters = [('A', 'C', 'G', 'T') if a == 'N' else (a,) for a in s]
    seqs = [''.join(s) for s in itertools.product(*letters)]
    return set(seqs)

def generate_similar(seq, m):
    """Generate sequences within some number of mismatches away from input.

    For a sequence in alphabet A, C, G, T, N, generate all sequences
    of the same length that have at most m mismatches with it, exluding any
    occurrences of N.

    Args:
    seq: str composed of A, C, G, T, N
    m: number of mismatches in generated sequences
    """
    # Note:
    # this function was written in hope that generating all these sequences
    # and then checking them in trie may be faster than trie.get_approximate
    # with the same mismatch parameter. But experiments showed it's not, it's
    # actually at least 2 times slower even just to check with trie.has_key(),
    # and on top of that this implementation of generate_similar is three times
    # slower than trie.get_approximate()
    comb = itertools.combinations(xrange(len(seq)), m)
    comb = ((('A', 'C', 'G', 'T') if i in x and c != 'N' else (c,)
             for i,c in enumerate(seq))
            for x in comb)
    comb = (itertools.product(*g) for g in comb)
    comb = (''.join(s) for s in itertools.chain(*comb))
    return comb

def map_coord_to_int(s, genome):
    """Map genomic coordinates to int.

    Validity of input is not checked. Chromosome name in 's' should be
    present in 'genome' and coordinate in 's' should be strictly less
    than length of this chromosome in 'genome' (0-based coordinate).

    Args:
    s: str of the form '<chromosome name>:<coordinate>:<strand>'
    genome: list of pairs [(<chromosome name>, <chromosome length>)]

    Return:
    int
    """
    chrom, coord, strand = s.split(':')
    coord = int(coord)
    x = 0
    i = 0
    while genome[i][0] != chrom:
        x += genome[i][1]
        i += 1
    x += coord
    x += 1  # to avoid problem with 0
    x = x if strand == '+' else -x
    return x

def map_int_to_coord(x, genome, onebased=False, strout=True):
    """Map int to genomic coordinates, inverse to map_coord_to_int().

    Validity of input is not checked. 'x' should be not greater than total
    length of all chromosomes in 'genome'. Output is 0-based coordinates
    by default unless parameter 'onebased' is True.

    Args:
    x: int to map
    genome: list of pairs [(<chromosome name>, <chromosome length>)]
    onebased: if True, return 1-based coordinates (standard for SAM etc.);
              otherwise stick to 0-based coordinates used internally
              in the trie etc.
    strout: if True, collapse output into str

    Return:
    tuple (chrom (str), coord (int), strand (str)) or
    str "<chromosome name>:<coordinate>:<strand>" depending on strout
    """
    strand = '+' if x > 0 else '-'
    x = abs(x)
    x -= 1  # see "x += 1" in reverse function
    i = 0
    while genome[i][1] <= x:
        x -= genome[i][1]
        i += 1
    chrom = genome[i][0]
    coord = x
    if onebased:
        coord += 1
    t = (chrom, coord, strand)
    if strout:
        return '%s:%s:%s' % t
    else:
        return t

def get_nonexist_int_coord(genome):
    """Calculate int coordinate that does not correspond to any real one.

    Args:
    genome: list of pairs [(<chromosome name>, <chromosome length>)]

    Return:
    int
    """
    return -(sum(p[1] for p in genome) + 1)

def array_to_hex(arr):
    """Map numpy.array of int to hex string.

    Note: may be system-dependent, need to test.
    """
    return binascii.hexlify(arr.tostring())

def hex_to_array(hexstr):
    """Map hex string to numpy.array of int.

    Note: may be system-dependent, need to test.
    """
    return np.fromstring(binascii.unhexlify(hexstr), dtype=int)

def offtargetinfo_to_hex(offtargetinfo, delim):
    """Map off-target information about guideRNA to hex string.

    Args:
    offtargetinfo: list of info about off-targets
                   [(similarity of off-target to guideRNA sequence,
                     array of a label and int coords for off-target)]
    delim: int that never appears in any of the arrays

    Return:
    hex str
    """
    arrays = [np.append(arr, [dist, delim])
              for dist, arr in offtargetinfo]
    return array_to_hex(np.concatenate(arrays))

def hex_to_offtargetinfo(hexstr, delim, offcoords=False):
    """Map hex string to off-target information about guideRNA.

    Args:
    hexstr: hex string that contains encoded off-target information
    delim: int that never appears in any of the actual arrays with int coords
           and labels
    offcoords: if True, return off-target coordinates (takes more time),
               otherwise return only distance and main label

    Return:
    shrinked off-target info in format
        [(distance, array of int coords and labels)]
    """
    mainarr = hex_to_array(hexstr)
    index = np.where(mainarr == delim)[0]
    index = np.concatenate(([-1], index))
    if offcoords:
        pairs = zip(index[:-1] + 1, index[1:])
        arrays = [(mainarr[p[1] - 1], mainarr[p[0]:p[1] - 1]) for p in pairs]
    else:
        arrays = zip(mainarr[index[1:] - 1], mainarr[index[:-1] + 1])
    return arrays
    # mainarr = hex_to_array(hexstr)
    # arrays = np.split(mainarr, np.where(mainarr == delim)[0] + 1)
    # arrays = [(arr[-2], arr[:-2]) for arr in arrays if len(arr) > 3]
    # return arrays

def create_inttree_from_file(infile):
    """Create interval tree to store annotations

    Args:
    infile: handle of open BED file with annotations

    Return:
    dictionary {chromosome name : interval tree with coordinates}
    """
    genome = {}
    for line in infile:
        clean_line = line.strip()
        parts = clean_line.split()
        chrom, start, stop = parts[0], int(parts[1]), int(parts[2])
        name = parts[3]
        tree = None
        #if chromosome already in tree, index to this tree
        if chrom in genome:
            tree = genome[chrom]
        else:
            #first time we encounter chromosome, create a new interval tree
            tree = IntervalTree()
            genome[chrom] = tree
        #add interval to tree
        tree.add(start, stop, name)
    return genome

def create_annot_inttree(annotfile):
    """Create interval tree with annotations.

    Args:
    annotfile: path to BED file with coordinates of genomic features
               in format (tabulated): chrom, start, end, name;
               or short name for preinstalled exon annotations
               in folder annotation_bed
    Return:
    dictionary {chromosome name : interval tree with coordinates}
    if successful or None otherwise
    """
    preinstalled_annotfile = "annotation_bed/%s/%s_exons.bed" \
                             % (annotfile, annotfile)
    if resource_exists(__name__, preinstalled_annotfile):
        inttree_filename = resource_filename(__name__, preinstalled_annotfile)
    else:
        inttree_filename = annotfile
    if not os.path.isfile(inttree_filename):
        return None
    infile = open(inttree_filename, 'r')
    genome = create_inttree_from_file(infile)
    infile.close()
    return genome
