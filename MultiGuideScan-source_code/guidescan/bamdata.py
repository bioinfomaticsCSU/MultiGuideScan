__author__ = 'Tao Li & Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Create and manipulate guideRNA database in BAM format.

produce_bam_custom(): main function, call subroutines in parallel
produce_bams_main(): default function when run from the main pipeline
trie_to_sam(): produce SAM file with guideRNA database from kmer trie
sam_to_bam(): produce sorted and indexed BAM from SAM
get_offtarget_info(): auxiliary function for trie_to_sam()
"""

#################
#               #
#   Libraries   #
#               #
#################

import argparse
import os, sys
from os.path import basename
import gzip
import tempfile
from hashlib import md5
from datetime import datetime

from Bio import Seq

import util
import guides
from multiprocessing import Pool,Process,Queue
#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def trie_to_sam(index, kmers_trie, keysfile, samfile, args, offdist, maxoffcount,
                process, n, parts):
    """Produce SAM file with guideRNAs and info about their off-targets.

    Convention: Off-target info is stored in the optional field with flag 'of'
    of type 'H' (hex byte array) produced by function
    util.offtargetinfo_to_hex() and can be restored using function
    util.hex_to_offtargetinfo(). Also store there optional fields 'od'
    and 'oc' of type 'i' (integer) indicating parameters offdist
    and maxoffcount, respectively, for which this off-target info was
    produced. This is needed if potentially different values of offdist
    and maxoffcount are used within the same SAM (i.e., for some
    guideRNAs one wants more detailed info about their off-targets
    than for others). 'od' and 'oc' are relevant even if 'of' is empty.
    If 'od' and 'oc' are not provided, refer to info from SAM header.


    Note: SAM uses 1-based genomic coordinates, and the SAM file produced
    by this function is consistent with that; 0-based coordinates of
    guideRNAs stored in the trie are transformed into 1-based coordinates.
    However, coordinates of off-targets stored in 'of' field remain
    0-based.

    Args:
    index: the index of keysfile
    kmers_trie: trie.trie object with all guideRNAs as produced by
                guides.analyze_guides()
    keysfile: name of file with all k-mers that are considered good
              candidate guideRNAs, one per line; if ends with '.gz'
              assume file is gzipped
    samfile: where to store SAM file, will be gzipped
    args: arguments of the project, used to print some info in SAM header
    offdist: maximum Hamming distance to consider from guideRNA to its
             off-target; use -1 for omitting any off-target info
             in resulting BAM (works much faster);
             use this value instead of what args contains
    maxoffcount: store at most this many off-targets for a guideRNA;
                 ignore if offdist is -1;
                 use this value instead of what args contains
    process: process number; to distinguish in output from different processes
    n: the number of prefix of kmers for preprocessing
    parts: the number of parts of classified kmers
    """
    util.print_log('process%s:trie_to_sam' % process)
    
    util.check_file_exists(keysfile)
    f = gzip.open(keysfile) if keysfile.endswith('.gz') else open(keysfile)
    s = gzip.open(samfile, 'w')
    s.write('@HD\tVN:1.0\tSO:unknown\n')
    genome = args['genome']
    # be careful with changing the next line
    # current rational is: arrays contain only global genomic coordinates
    # (and label as 0-th element); delim avoids all of these
    delim = util.get_nonexist_int_coord(genome)
    for chrom, length in genome:
        s.write('@SQ\tSN:%s\tLN:%s\n' % (chrom, length))
    s.write('@CO\tprepared with iGuide software\n')
    s.write('@CO\tcontains info about guideRNAs and their off-targets\n')
    s.write('@CO\targuments of the run:\n')
    s.write('@CO\t%s\n' % args)
    count = 0
    starttime = datetime.now()
    lasttime = starttime

    index_seq = guides.generate_four(index, n)
    mismatches = []
    
    if offdist != -1:
        for i in range(parts):
            index_seq1 = guides.generate_four(i, n)
            mismatch = guides.four_compare(index_seq, index_seq1, n)
            mismatches.append(mismatch)
    # util.print_log('process%s:mismatch done...' % process)
    for line in f:
        
        guide = line.split()[0]
        kmer2 = guide[n:]
        
        count += 1
        
        # QNAME
        samline = guide
#        samline = '%st%s' % (count, process)
        if not kmers_trie[index].has_key(kmer2):
            print 'process %s warning: %s is not in trie %s, skip' \
                  % (process, kmer2, index)
            continue
        arr = kmers_trie[index][kmer2]
        if arr[0] != 0:
            print 'process %s warning: %s is not a good guideRNA according' \
                  ' to label in trie %s, skip' % (process, kmer2, index)
            continue
        if len(arr) > 2:
            print 'process %s warning: %s is stored with more than one' \
                  ' coordinate in trie %s, skip' % (process, kmer2, index)
            continue
        coord = arr[1]
        coord = util.map_int_to_coord(coord, genome)
        chrom, pos, strand = coord.split(':')
        pos = int(pos)
        flag = '0' if strand == '+' else '16'
        # FLAG
        samline += '\t%s' % flag
        # RNAME
        samline += '\t%s' % chrom
        if strand == '-':
            pos = pos - args['length'] - len(args['pam']) + 1
        pos += 1  # SAM uses 1-based coordinates, in our code we used 0-based
        # POS
        samline += '\t%s' % pos
        # MAPQ
        samline += '\t100'  # 100 is arbitrary choice, 255 is not recommended
        # CIGAR
        samline += '\t%sM' % (args['length'] + len(args['pam']))
        # RNEXT
        samline += '\t*'
        # PNEXT
        samline += '\t0'
        # TLEN
        samline += '\t0'
        # SEQ
        seq = guide
        if strand == '-':
            seq = str(Seq.Seq(guide).reverse_complement())
        samline += '\t%s' % seq
        # QUAL
        samline += '\t*'
        # offtargets
        offtargetargs = 'od:i:%s\toc:i:%s' % (offdist, maxoffcount)
        samline += '\t%s' % offtargetargs
        if offdist != -1:

            offtargetinfo = get_offtarget_info(mismatches, kmers_trie, kmer2, offdist,
                                               maxoffcount, delim, parts)
            if offtargetinfo:
                offtargetstr = 'of:H:%s' % offtargetinfo
                samline += '\t%s' % offtargetstr
        samline += '\n'
        s.write(samline)
        currenttime = datetime.now()
        # print every minute during first 10 minutes, then every 20 minutes
        if (currenttime - lasttime).seconds > 1200 \
           or ((currenttime - lasttime).seconds > 60
               and (lasttime - starttime).seconds < 600):
            util.print_log('process %s: %s guides processed' % (process, count))
            lasttime = currenttime
    util.print_log('process %s: total %s guides processed' % (process, count))
    s.flush()
    s.close()
    f.close()

def get_offtarget_info(mismatches, kmers_trie, kmer2, offdist, maxoffcount, delim, parts):
    """Get off-target info about a guideRNA and output as str.

    Args:
    kmers_trie: trie.trie object with all guideRNAs as produced by
                guides.analyze_guides()
    guide: sequence of guideRNA
    offdist: maximum Hamming distance to consider from guideRNA to
             its off-target; running time icreases somewhat exponentially
             as this value increases; offdist=4 may be infeasible
             when running genome-wide analysis on mammalian genome
    maxoffcount: store at most this many off-targets for a guideRNA
    delim: int that never appears in any of the arrays at values of kmers_trie

    Return:
    hex str with off-target info (as produced by util.offtargetinfo_to_hex())
    """
    # conclusion after some testing: the next call
    # to trie.get_approximate_hamming() is the bottleneck of the whole
    # SAM file production;
    # offdist=3 seems feasible but larger values may be infeasible when
    # looping over millions of guideRNAs and trie built from really
    # large genome
    offtargets = []
    for i in range(parts):
        mismatch = mismatches[i]
        if mismatch > offdist:
            continue
        offdist1 = offdist - mismatch
        offtarget1 = kmers_trie[i].get_approximate_hamming(kmer2, offdist1)
        if (len(offtarget1) < 1):
            continue
        offtarget2 = [(offtarget, info, distance+mismatch)
                        for offtarget, info, distance in offtarget1]
        offtargets.extend(offtarget2)

    if len(offtargets) <= 1:
        return
    # in next line, sort by distance, then by count of occurrences,
    # and then within the same distance reproducibly shuffle at random
    # to avoid bias in reporting
    offtargets = sorted(offtargets,
                        key = lambda p: (p[2], -p[1][0], md5(p[0]).digest()))
    offtargetinfo = [(distance, info)
                     for offtarget, info, distance in offtargets
                     if distance > 0]
    if not offtargetinfo:
        return
    offtargetinfo = offtargetinfo[:maxoffcount]
    offtargetinfo = util.offtargetinfo_to_hex(offtargetinfo, delim)
    return offtargetinfo

def sam_to_bam(samfile, bamfile, index=False):
    """Produce sorted and indexed BAM file from SAM file with guideRNAs.

    Relies on 'gzip', 'samtools' available in the system.

    Args:
    samfile: where SAM file is stored, assume file is gzipped
    bamfile: where to store BAM file
    index: if True, index the resulting BAM file
    """
    util.check_file_exists(samfile)
    # util.warn_file_exists(bamfile)
    samtools_command = 'gzip -cd %s | samtools view -hb - ' \
                       '| samtools sort - > %s' \
                       % (samfile, bamfile)
    # print samtools_command
    os.system(samtools_command)
    if index:
        samtools_index_command = 'samtools index %s' % bamfile
        print samtools_index_command
        os.system(samtools_index_command)

def process_pool(q, kmers_trie_list, args, offdist, maxoffcount, process, n, parts):

    while True:
        try:
            param = q.get(False)
            keysfile = param[0]
            samfile = param[1]
            index = param[2]
            
            trie_to_sam(index, kmers_trie_list, keysfile, samfile, args, offdist, maxoffcount, process, n, parts)
        except Exception:
            if q.empty():
                break

def produce_bam_custom(kmers_trie, name, label, guides_filename, args,
                       offdist, maxoffcount, processes, n, parts):
    """Produce BAM file with guideRNA database.

    Run after all files and trie were generated
    by kmers.extract_process_kmers() and guides.analyze_guides()

    Produce files:
    sorted BAM file with off-target info:
        <name>/<name>_guides_<label>.bam
    index for the BAM file with off-target info:
        <name>/<name>_guides_<label>.bam.bai

    Args:
    kmers_trie: trie.trie object with all guideRNAs as produced by
                guides.analyze_guides()
    name: project name, used to get project args and in all output
    label: str, add it to file name of output database for this run
    guides_filename: name of file with all k-mers that are considered good
                     candidate guideRNAs, one per line;
                     if file name ends with .gz assume file is gzipped;
    args: arguments of the project, used to print some info in SAM header
    offdist: maximum Hamming distance to consider from guideRNA to
             its off-target;
             use -1 for omitting any off-target info in resulting BAM
             (works much faster);
             running time icreases somewhat exponentially
             as this value increases; offdist=4 may be infeasible
             when running genome-wide analysis on mammalian genome
    maxoffcount: store at most this many off-targets for a guideRNA;
                 ignore if offdist is -1
    processes: int, how many processes to use in parallel; do not specify more
             than available in the system; currently not implemented, use 1
    """
    guidesfiles = []
    # parts = 256
    tempdir = '%s%s' % (name,'/classifiedfiles/tempfiles')

    util.print_log('produce SAM files...')
    samfiles = ['%s/%s.sam' % (tempdir, i) for i in range(parts)]
    # samfiles = [tempfile.NamedTemporaryFile(dir=name, suffix='.sam%s' % i)
    #             for i in xrange(parts)]
    # util.print_log('store SAM in these files (gzipped): %s'
    #                % (', '.join([basename(f.name) for f in samfiles])))

        
    if isinstance(guides_filename, str):

        util.print_log('split %s in %s parts...' % (guides_filename, parts))
        guidesfiles = [tempfile.NamedTemporaryFile(dir=name,
                                                   suffix='.guides%s' % i)
                       for i in range(parts)]
        util.print_log('store guides in these files: %s'
                       % (', '.join([basename(f.name) for f in guidesfiles])))
        guidesfile = gzip.open(guides_filename) \
                     if guides_filename.endswith('.gz') \
                     else open(guides_filename)
        index_num = 0
        guidecount = 0
        for line in guidesfile:
            kmer1 = line.split()[0][0:n]
            index_num = guides.get_num(kmer1, n)
            guidesfiles[index_num].write(line)
            
            guidecount += 1
            
        guidesfile.close()
        for f in guidesfiles:
            f.flush()
        util.print_log('%s guideRNAs to process' % guidecount)
        util.print_log('done')

        process_list = []
        all_task = Queue()
        for i in range(parts):
            task = (guides_filename[i].name, samfiles[i].name, i)
            all_task.put(task)

        for i in range(processes):
            p = Process(target=process_pool, args=(all_task, kmers_trie, args, offdist, maxoffcount, i, n, parts))
            p.start()
            process_list.append(p)

        for p in process_list:
            p.join()
                
        for i in range(parts):
            guidesfiles[i].close()

    else:       
        process_list = []
        all_task = Queue()
        for i in range(parts):
            task = (guides_filename[i], samfiles[i], i)
            all_task.put(task)

        for i in range(processes):
            p = Process(target=process_pool, args=(all_task, kmers_trie, args, offdist, maxoffcount, i, n, parts))
            p.start()
            process_list.append(p)

        for p in process_list:
            p.join()

    util.print_log('produce sorted BAM files...')
    
    bamfiles = ['%s/%s.bam' % (tempdir, i) for i in range(parts)]
    # bamfiles = [tempfile.NamedTemporaryFile(dir=name, suffix='.bam%s' % i)
    #             for i in xrange(parts)]
    # util.print_log('store BAM in these files: %s'
    #                % (', '.join([basename(f.name) for f in bamfiles])))

    pool = Pool(processes)
    util.print_log('poolSize %s...' % processes)
    index=False
    for i in range(parts):
        pool.apply_async(sam_to_bam,(samfiles[i], bamfiles[i], index,))
    util.print_log('Waiting for all subprocesses done...')
    pool.close()
    pool.join()

    # for i in xrange(parts):
    #     samfiles[i].close()
    util.print_log('merge into one BAM file...')
    bamfile = '%s/%s_guides%s.bam' % (name, name,
                                      '_%s' % label if label else '')
    util.print_log('store in %s' % bamfile)
    util.warn_file_exists(bamfile)
    if parts > 1000:
        mid = parts // 2
        bamfiles_temp = [tempfile.NamedTemporaryFile(dir=name, suffix='.bam%s' % i)
                        for i in xrange(2)]
        samtools_command1 = 'samtools merge -f %s %s' \
                           % (bamfiles_temp[0].name, ' '.join(bamfiles[0:mid]))
        os.system(samtools_command1)

        samtools_command2 = 'samtools merge -f %s %s' \
                           % (bamfiles_temp[1].name, ' '.join(bamfiles[mid:parts]))
        os.system(samtools_command2)

        samtools_command = 'samtools merge -f %s %s' \
                           % (bamfile, ' '.join([f.name for f in bamfiles_temp]))
        os.system(samtools_command)

        for f in bamfiles_temp:
            f.close()

    else:
        samtools_command = 'samtools merge -f %s %s' \
                           % (bamfile, ' '.join(bamfiles))
        # print samtools_command
        os.system(samtools_command)
    samtools_index_command = 'samtools index %s' % bamfile
    # print samtools_index_command
    os.system(samtools_index_command)
    util.print_log('done')
    # for i in xrange(parts):
    #     bamfiles[i].close()

    for i in range(parts):
        if(os.path.exists(samfiles[i])):
            os.remove(samfiles[i])
        if(os.path.exists(bamfiles[i])):
            os.remove(bamfiles[i])

    util.print_log('samtools version')
    samtools_version_command = 'samtools --version'
    print samtools_version_command
    os.system(samtools_version_command)

def produce_bams_main(kmers_trie, name):
    """Produce BAM file with all guideRNAs and info about their off-targets.

    Run after all files and trie were generated
    by kmers.extract_process_kmers() and guides.analyze_guides()

    Produce files:
    sorted BAM file with off-target info: <name>/<name>_guides.bam
    index for the BAM file with off-target info: <name>/<name>_guides.bam.bai
    also, BAM file and index for all guideRNAs without any off-target info
    (produced much faster):
        <name>/<name>_guides_nooff.bam
        <name>/<name>_guides_nooff.bam.bai

    Args:
    kmers_trie: trie.trie object as produced by guides.analyze_guides()
    name: project name, used to get project args and in all output
    """
    util.print_log('start produce_bam()')
    util.print_log('load arguments...')
    args = util.load_args(name)
    util.print_args(args)
    util.print_log('done')

    util.print_log('produce SAM file with guideRNAs only (no off-targets)...')
    # guides_filename = '%s/%s_guides.txt.gz' % (name, name)
    # parts = 256
    n = args['greateroffdist']
    parts = 4 ** n

    guides_dir = '%s%s' % (name,'/classifiedfiles/guides')
    guides_filenames = ['%s/%s.txt.gz' % (guides_dir, i) for i in range(parts)]
    util.print_log('read guides from %s' % guides_dir)
    produce_bam_custom(kmers_trie=kmers_trie, name=name, label='nooff',
                       guides_filename=guides_filenames,
                       args=args, offdist=-1,  # -1 for no off-targets
                       maxoffcount=args['maxoffcount'],
                       processes=args['processes'],
                       n = n,
                       parts=parts)
    util.print_log('done')

    if args['offdist'] != -1:
        util.print_log('produce SAM file with guideRNAs'
                       ' and off-target info...')
        # guides_filename = '%s/%s_guides.txt.gz' % (name, name)
        util.print_log('read guides from %s' % guides_dir)
        produce_bam_custom(kmers_trie=kmers_trie, name=name, 
                           label='offdist%s' % args['offdist'],
                           guides_filename=guides_filenames,
                           args=args, offdist=args['offdist'],
                           maxoffcount=args['maxoffcount'],
                           processes=args['processes'],
                           n = n,
                           parts=parts)
        util.print_log('done')

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():
    p = argparse.ArgumentParser(description='Produce BAM file with guideRNA'
                                            ' database from precomputed trie'
                                            ' and list of guideRNAs',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-n', dest='name', default='myguides', required=True,
                   help='project name, load previously saved arguments'
                        ' and save additional output')
    p.add_argument('--label', dest='label', default='test', required=True,
                   help='use in file name of output database for this run')
    p.add_argument('-g', dest='guidesfile', default='',
                   help='name of file with guideRNAs for which to compute'
                        ' BAM database; may be gzipped (.gz);'
                        ' if not provided, use all candidate guideRNAs'
                        ' found in the project')
    p.add_argument('-d', dest='offdist', type=int, default=3,
                   help='maximum Hamming distance to consider from guideRNA'
                        ' to its off-target;'
                        ' off-target is an alternative occurrence (with any'
                        ' PAM) of this guideRNA in the genome at Hamming' 
                        ' distance at most this number (including PAM);'
                        ' use -1 for omitting any off-target info in resulting'
                        ' BAM (works much faster)')
    p.add_argument('-k', dest='greateroffdist', type=int, default=4,
                   help='a number greater than offdist used for preprocessed data'
                        '(the length of key for classifying guide RNAs)')
    p.add_argument('--maxoffcount', dest='maxoffcount', type=int, default=1000,
                   help='maximum number of off-targets to store for'
                        ' a guideRNA in a resulting BAM library;'
                        ' ignore if OFFDIST is -1')
    p.add_argument('-t', dest='processes', type=int, default=1,
                   help='how many processes to use; do not specify more'
                        ' than you have on your system;'
                        ' currently not implemented')
    args = p.parse_args()
    sam_args_dict = args.__dict__
    name = sam_args_dict['name']
    guides_filename = sam_args_dict['guidesfile']

    n = sam_args_dict['greateroffdist']
    parts = 4 ** n
    # parts = 256
    if not guides_filename:
        # guides_filename = '%s/%s_guides.txt.gz' % (name, name)
        guides_dir = '%s%s' % (name,'/classifiedfiles/guides')
        guides_filename = ['%s/%s.txt.gz' % (guides_dir, i) for i in range(parts)]

    util.print_log('local script arguments:')
    util.print_args(sam_args_dict)
    util.print_log('load main arguments...')
    args = util.load_args(name)
    util.print_args(args)
    util.print_log('done')

    # main
    trie_filename = ['%s/%s/%s_trie%s.dat' % (name, 'kmers_tries', name, i) for i in range(parts)]
    kmers_trie = guides.load_restore_trie(name, trie_filename, n, parts)
    produce_bam_custom(kmers_trie=kmers_trie,
                       name=name,
                       label=sam_args_dict['label'],
                       guides_filename=guides_filename,
                       args=args,
                       offdist=sam_args_dict['offdist'],
                       maxoffcount=sam_args_dict['maxoffcount'],
                       processes=sam_args_dict['processes'],
                       n = n,
                       parts=parts)


if __name__ == '__main__':
    main()
