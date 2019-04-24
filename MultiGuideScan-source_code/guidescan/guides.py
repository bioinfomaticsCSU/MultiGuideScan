# -*- coding: UTF-8 -*-
__author__ = 'Alexendar Perez & Yuri Pritykin'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""Analyze k-mers and extract guideRNAs with all relevant info about them.

analyze_guides(): main function, create, modify and save trie, print log
build_kmers_trie(): create initial trie of k-mers and their coordinates
label_multimapping(): label multimapping k-mers in trie with correct counts
filter_trie_mismatch_similarity(): find keys in trie similar to other keys
filter_keys_trie(): using labels in trie, choose good k-mers from a list

trie_get_approximate(): wrapper of trie.trie.get_approximate()
trie_get_approximate_hamming(): deprecated; use trie.get_approximate_hamming()

save_trie(): save trie to file
load_trie(): load trie from file
restore_trie_arrays(): restore numpy.array from str in loaded kmer trie
load_restore_trie(): load and restore trie of all genomic kmers for a project
"""

#################
#               #
#   Libraries   #
#               #
#################

import os,sys
import argparse
from os.path import basename
import gzip
from collections import defaultdict
import tempfile
import numpy as np
# from Bio import trie

import util
import trie
from multiprocessing import Pool,Process,Queue,Lock
import multiprocessing
#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

#Converting quaternary string to decimal
def get_num(kmer, n):
    four = kmer[0:n]
    ret = 0
    for c in four:
        ret *= 4
        if c == 'A': ret += 0
        if c == 'C': ret += 1
        if c == 'G': ret += 2
        if c == 'T': ret += 3
    return ret
#Converting decimal number to quaternary string
def generate_four(num, n):
    ret = ''
    for i in range(0, n):
        p = num % 4
        if p == 0: ret = 'A' + ret
        if p == 1: ret = 'C' + ret
        if p == 2: ret = 'G' + ret
        if p == 3: ret = 'T' + ret
        num //= 4
    return ret

def four_compare(s1, s2, n):
    count = 0
    for i in range(0, n):
        if s1[i] != s2[i]:
            count += 1
    return count

def save_trie(t, filename, parts):
    """Save trie to file.

    Note: saving and loading tries (especially with np.array values)
    may be system-dependent (e.g., depend on 64bit or 32bit arithmetic),
    need to test.

    Args:
    t: trie.trie object
    filename: where to save
    """
    for i in range(parts):
        f = open(filename[i], 'w')
        trie.save(f, t[i])
        f.close()

def save_single_trie(t, filename):
    """Save trie to file.

    Note: saving and loading tries (especially with np.array values)
    may be system-dependent (e.g., depend on 64bit or 32bit arithmetic),
    need to test.

    Args:
    t: trie.trie object
    filename: where to save
    """
    
    f = open(filename, 'w')
    trie.save(f, t)
    f.close()

def load_trie(filename, parts):
    """Load trie from file.

    Note: saving and loading tries (especially with np.array values)
    may be system-dependent (e.g., depend on 64bit or 32bit arithmetic),
    need to test.
    """
    kmers_trie = []
    for i in range(parts):
        f = open(filename[i])
        t = trie.load(f)
        kmers_trie.append(t)
        f.close()
    return kmers_trie

def restore_trie_arrays(kmers_trie, keysfile, n):
    """Restore numpy arrays from str in values of kmer trie loaded from disk.

    When loaded from disk, numpy arrays in values of a previously stored trie
    are loaded as string characters. This function restores the initial arrays
    using numpy.fromstring(). To save memory, avoid building a list of keys
    of the trie and take them from file.

    Note: saving and loading tries (especially with np.array values)
    may be system-dependent (e.g., depend on 64bit or 32bit arithmetic
    and endianness), use caution when transferring files between systems.

    Args:
    kmers_trie: trie.trie object loaded from dist (e.g., using load_trie())
    keysfile: name of file where first field of each line is a k-mer, assume
              file is gzipped; loop only over the k-mers in this file
    """
    util.check_file_exists(keysfile)
    f = gzip.open(keysfile)
    for line in f:
        kmer = line.split()[0]
        index = get_num(kmer, n)
        kmer2 = kmer[n:]

        if kmers_trie[index].has_key(kmer2):
            value = kmers_trie[index][kmer2]
            if isinstance(value, basestring):
                kmers_trie[index][kmer2] = np.fromstring(value, dtype=int)
    f.close()
    return kmers_trie

def load_restore_trie(name, trie_filename, n, parts):
    """Fully load previously stored trie of all genomic kmers for a project.

    name: project name, used to get project args and in all output

    Return:
    loaded kmers_trie
    """
    util.print_log('load trie...')
    # trie_filename = ['%s/%s/%s_trie%s.dat' % (name, 'kmers_tries', name, i) for i in range(256)]
    kmers_trie = load_trie(trie_filename, parts)
    util.print_log('done')
    keysfile = '%s/%s_kmers_shuffled.txt.gz' % (name, name)
    util.print_log('restore numpy arrays in trie values...')
    restore_trie_arrays(kmers_trie, keysfile, n)
    util.print_log('done')
    return kmers_trie

def build_kmers_trie(filename, genome, name, altpam=[], pampos='end', maxcount=10,
                     goodkeysfile='', badkeysfile='', tempdir='', triekeys_v1_filenames=[],
                      kmers_filenames=[], processes=1, n=4, parts=256):
    """Read k-mers and their coordinates from file and store them in a trie.

    The resulting trie is of the form {<k-mer> : <np.array of int values>}.
    In each array, store int-transformed coordinates of occurrences of
    the k-mer, starting from position 1 of the array. Store at most
    'maxcount' coordinates.
    Position 0 in the array is reserved for labeling:
    0: good guideRNA,
    positive: show how many occurrences this k-mer has in the genome
    [other labels: decide later]
    In this building stage, label all k-mers with alternative PAM
    and all k-mers with more than one occurrence in the genome
    as bad guideRNAs.
    Optionally, also store in a separate file all k-mers that are still
    considered good candidate guideRNAs. This means only filtering out
    k-mers with alternative PAM, because detecting multi-mapping k-mers
    labeling them as bad guideRNAs may happen after they were first read.

    Note: make sure lines with k-mers in input file are randomly shuffled.
    This is to ensure that for k-mers with more than 'maxcount' occurrences,
    we store artibrary 'maxcount' of them without any bias.

    Args:
    filename: name of file with k-mers, assume file is gzipped and
              lines are randomly shuffled
    genome: list of pairs [(<chromosome name>, <chromosome length>)]
    altpam: list of alternative PAM sequences, all k-mers starting
            or ending (depending on argument 'pampos') with these
            sequences are labeled as bad guideRNAs
    pampos: position of alternative PAM in k-mer ('start' or 'end')
    maxcount: store at most this many coordinates for each k-mer
    goodkeysfile: where to store potentially good candidate guideRNAs;
                  use only if altpam is not empty, otherwise all input
                  keys from filename will be stored which is redundant

    Output:
    return trie.trie object {<k-mer> : <np.array of int values>}
    optionally produce file goodkeysfile with candidate guideRNAs
    """
    # parts = 256
    util.check_file_exists(filename)   

    badkeysfiles = ['%s/badkeys%s.txt.gz' % (tempdir, i) for i in range(parts)]

    kmers_trie_files = ['%s/kmers_trie%s.dat' % (tempdir, i) for i in range(parts)]

    util.print_log('classify k-mers into %s...' % parts)
    if parts > 1000:
        tempfiles = [tempfile.NamedTemporaryFile(dir=tempdir,
                                               suffix='.temp%s' % i)
                        for i in range(2)]
        # tempfiles = [gzip.open('%s/temp%s.txt.gz' % (tempdir, i),'w') for i in range(2)]
        mid = parts // 2
        file = gzip.open(filename)
        for line in file:
            kmer = line.split()[0]
            index = get_num(kmer, n)
            if index < mid:
                tempfiles[0].write(str(index) + ' ' + line)
            else:
                tempfiles[1].write(str(index) + ' ' + line)
        for f in tempfiles:
            f.flush()
        file.close()
        util.print_log('write...')

        kmersfiles1 = [gzip.open(kmers_filenames[i],'w') for i in range(mid)]
        # tempfiles = [gzip.open('%s/temp%s.txt.gz' % (tempdir, i)) for i in range(2)]
        temp = [open(tempfiles[i].name) for i in range(2)]
        for line in temp[0]:
            data = line.split()
            index = int(data[0])
            kmer = data[1]
            coord = data[2]
            kmersfiles1[index].write(kmer + ' ' + coord + '\n')
        for f in kmersfiles1:
            f.close()
        temp[0].close()
        util.print_log('write count1...')

        kmersfiles2 = [gzip.open(kmers_filenames[i],'w') for i in range(mid, parts)]
        for line in temp[1]:
            data = line.split()
            index = int(data[0]) - mid
            kmer = data[1]
            coord = data[2]
            kmersfiles2[index].write(kmer + ' ' + coord + '\n')
        for f in kmersfiles2:
            f.close()
        temp[1].close()
        for f in tempfiles:
            f.close()       
        util.print_log('write count2...')
    else:
        kmersfiles = [gzip.open(kmers_filenames[i],'w') for i in range(parts)]
        file = gzip.open(filename)
        
        for line in file:
            kmer = line.split()[0]
            index = get_num(kmer, n)
            kmersfiles[index].write(line) 
        file.close()
        for f in kmersfiles:
            f.close()

    util.print_log('done...')
    util.print_log('build tries start...')
    
    process_list = []
    all_task = Queue()
    for i in range(parts):
        task = (kmers_filenames[i], triekeys_v1_filenames[i], badkeysfiles[i], kmers_trie_files[i])
        all_task.put(task)

    for process in range(processes):
        p = Process(target=process_pool_build_tries, args=(all_task, genome, altpam, pampos, maxcount, n))
        p.start()
        process_list.append(p)

    for p in process_list:
        p.join()

    util.print_log('build tries done...')

    if goodkeysfile:
        util.warn_file_exists(goodkeysfile)
        goodkeys = gzip.open(goodkeysfile, 'w')
    if badkeysfile:
        util.warn_file_exists(badkeysfile)
        badkeys = gzip.open(badkeysfile,'w')

    count_added = 0
    for i in range(parts):
        keys_filename = triekeys_v1_filenames[i]
        f = gzip.open(keys_filename)
        for line in f:
            goodkeys.write(line)
            count_added += 1
        f.close()
    for i in range(parts):
        f = gzip.open(badkeysfiles[i])
        for line in f:
            badkeys.write(line)
            count_added += 1
        f.close()
    goodkeys.close()
    badkeys.close()
    
    # kmers_trie_filenames = [kmers_trie_files[i].name for i in range(parts)]
    kmers_trie = load_restore_trie(name, kmers_trie_files, n, parts)
    for i in range(parts):
        if(os.path.exists(kmers_trie_files[i])):
            os.remove(kmers_trie_files[i])
    print '%s keys added to trie' % count_added
    return kmers_trie

def process_pool_build_tries(queue, genome, altpam, pampos, maxcount, n):

    while True:
        try:
            param = queue.get(False)
            kmers_filename = param[0]
            goodkeys_filename = param[1]
            badkeys_filename = param[2]
            kmers_trie_filename = param[3]
            
            build_kmers_tries(kmers_filename, goodkeys_filename, badkeys_filename, 
                              kmers_trie_filename, genome, altpam, pampos, maxcount, n)
        except Exception, e:
            # print "Exception occured: %s" % e
            if queue.empty():
                break

def build_kmers_tries(kmers_filename, goodkeys_filename, badkeys_filename, 
                      kmers_trie_filename, genome, altpam, pampos, maxcount, n):
    util.check_file_exists(kmers_filename)
    if goodkeys_filename:
        goodkeys = gzip.open(goodkeys_filename, 'w')
    if badkeys_filename:
        badkeys = gzip.open(badkeys_filename,'w')

    kmers_trie = trie.trie()

    f = gzip.open(kmers_filename)
    for line in f:
        kmer, coord = line.strip().split()
        kmer2 = kmer[n:]

        if kmers_trie.has_key(kmer2):
            arr = kmers_trie[kmer2]
            if len(arr) < maxcount + 1:
                coord_int = util.map_coord_to_int(coord, genome)
                arr = np.append(arr, coord_int)
                arr[0] = len(arr) - 1
                kmers_trie[kmer2] = arr
        else:
            coord_int = util.map_coord_to_int(coord, genome)
            label = 0
            if pampos == 'start' and any(kmer.startswith(p) for p in altpam):
                label = 1
            if pampos == 'end' and any(kmer.endswith(p) for p in altpam):
                label = 1
            kmers_trie[kmer2] = np.array([label, coord_int])
            # count_added.value += 1
            if label == 0:
                goodkeys.write('%s\n' % kmer)
            if label != 0:
                badkeys.write('%s\n' % kmer)
    
    save_single_trie(kmers_trie, kmers_trie_filename)

    goodkeys.close()
    badkeys.close()
    f.close()

def label_multimapping(kmers_trie, filename, n):
    """Read multimapping k-mers and counts from file and add counts to trie.

    For each kmer and its count, find it in input trie and update the count
    of the kmer stored in position 0 of the array in the value of the trie.

    Args:
    kmers_trie: trie.trie object of the form described in and returned
                by build_kmers_trie()
    filename: name of file with k-mers and counts for k-mers occurring too
              many times in the genome,
              e.g., file <name>/<name>_kmers_counts.txt.gz
              produced by kmers.sort_count_kmers();
              assume file is gzipped
    
    Return:
    modified input trie.trie object
    """
    util.check_file_exists(filename)
    f = gzip.open(filename)
    count_labeled = 0
    for line in f:
        kmer, count = line.split()
        index = get_num(kmer, n)
        kmer2 = kmer[n:]

        if kmers_trie[index].has_key(kmer2):
            kmers_trie[index][kmer2][0] = count
            count_labeled += 1
    f.close()
    print '%s k-mers assigned counts' % count_labeled
    return kmers_trie

def filter_keys_trie(tempdir, kmers_trie, filenames1, filenames2, keysoutputfile,
                     nonCandidatekeysoutputfile, processes, n, parts):
    """Select k-mers from file that have label 0 in trie and write to file.

    Args:
    kmers_trie: trie.trie object of the form described in and returned
                by build_kmers_trie()
    keysinputfile: name of file where first field of each line is a k-mer,
                   assume file is gzipped
    keysoutputfile: name of file where to write selected keys one per line,
                    gzipped
    """
    # parts = 256
    badkeysfiles = ['%s/badkeys%s.txt.gz' % (tempdir, i) for i in range(parts)]
    
    process_list = []
    all_task = Queue()
    for i in range(parts):
        task = (filenames1[i], filenames2[i], badkeysfiles[i], i)
        all_task.put(task)

    for process in range(processes):
        p = Process(target=process_pool_filter, args=(all_task, kmers_trie, n))
        p.start()
        process_list.append(p)

    for p in process_list:
        p.join()  
    
    util.print_log('filter processes done...')
    badkeys = gzip.open(nonCandidatekeysoutputfile, 'w')
    for i in range(parts):
        f = gzip.open(badkeysfiles[i])
        for line in badkeysfiles[i]:
            badkeys.write(line)
        f.close()
    badkeys.close()

    write_count = 0
    goodkeys = gzip.open(keysoutputfile, 'w')
    for i in range(parts):
        f = gzip.open(filenames2[i])
        for line in f:
            goodkeys.write(line)
            write_count += 1
        f.close()

    print '%s keys written' % write_count

def process_pool_filter(queue, kmers_trie_list, n):

    while True:
        try:
            param = queue.get(False)
            keysinputfile = param[0]
            keysoutputfile = param[1]
            badkeysfile = param[2]
            index = param[3]
            
            filter_keys(kmers_trie_list, keysinputfile, keysoutputfile, badkeysfile, index, n)
        except Exception, e:
            # print "Exception occured: %s" % e
            if queue.empty():
                break

def filter_keys(kmers_trie, keysinputfile, keysoutputfile, badkeysfile, index, n):
    
    keysinput = gzip.open(keysinputfile)   
    keysoutput = gzip.open(keysoutputfile, 'w')
    badkeys = gzip.open(badkeysfile, 'w')
    for line in keysinput:
        kmer = line.split()[0]
        kmer2 = kmer[n:]

        if kmers_trie[index].has_key(kmer2) and kmers_trie[index][kmer2][0] == 0:
            keysoutput.write('%s\n' % kmer)
        elif kmers_trie[index].has_key(kmer2) and kmers_trie[index][kmer2][0] != 0:
            badkeys.write('%s\n' % kmer)
    keysoutput.close()
    badkeys.close()
    keysinput.close()


   # store k-mers whose label should be modified and finally modify tries together
def filter_trie_mismatch_similarity(tempdir, name, kmers_trie, sim, keys_filenames, processes, n, parts):
    """Find keys in trie few mismatches away from other keys.

    Args:
    kmers_trie: trie.trie object of the form described in and returned
                by build_kmers_trie()
    sim: if a key has another key at Hamming distance at most this,
         label it as bad guideRNA
    keysfile: name of file where first field of each line is a k-mer, assume
              file is gzipped; loop only over the k-mers in this file
    """
    
    # parts = 256
    badkeysfiles = ['%s/badkeys%s.txt.gz' % (tempdir, i) for i in range(parts)]
    
    process_list = []
    all_task = Queue()
    for i in range(parts):
        task = (keys_filenames[i], i, badkeysfiles[i])
        all_task.put(task)

    for process in range(processes):
        p = Process(target=process_pool_mismatch, args=(all_task, kmers_trie, sim, n, parts))
        p.start()
        process_list.append(p)

    for p in process_list:
        p.join()
 
    count = 0
    for i in range(parts):
        badkeys = gzip.open(badkeysfiles[i])

        for line in badkeys:
            string = line.strip().split()
            index = int(string[0])
            kmer = string[1]
            if not kmers_trie[index].has_key(kmer):
                continue
            if kmers_trie[index][kmer][0] != 0:
                continue
            else:
                kmers_trie[index][kmer][0] = 1
                count += 1
        badkeys.close()
    print '%s k-mers labeled as bad guideRNAs' % count
    util.print_log('done')
    
    return kmers_trie

def process_pool_mismatch(queue, kmers_trie_list, sim, n, parts):

    while True:
        try:
            param = queue.get(False)
            keysfile = param[0]
            index = param[1]
            badkeysfile = param[2]
            
            filter_trie_mismatch(badkeysfile, index, kmers_trie_list, keysfile, sim, n, parts)
        except Exception, e:
            # print "Exception occured: %s" % e
            if queue.empty():
                break

def filter_trie_mismatch(badkeysfile, index, kmers_trie, keysfile, sim, n, parts):
    
    util.check_file_exists(keysfile)

    index_seq = generate_four(index, n)
    mismatches = []
    for i in range(parts):
        index_seq1 = generate_four(i, n)
        mismatch = four_compare(index_seq, index_seq1, n)
        mismatches.append(mismatch)

    f = gzip.open(keysfile)

    badkeys = gzip.open(badkeysfile, 'w')
    
    for line in f:
        kmer = line.split()[0]
        kmer2 = kmer[n:]

        if not kmers_trie[index].has_key(kmer2):
            continue
        if kmers_trie[index][kmer2][0] != 0:
            continue 

        for i in range(parts):
            
            mismatch = mismatches[i]
            if mismatch > sim:
                continue 
            else:
                all_dist = kmers_trie[i].get_approximate_hamming(kmer2, sim-mismatch)

                if any(dist+mismatch > 0 for seq,arr,dist in all_dist):
                    badkeys.write('%s %s\n' % (index, kmer2))
                
    f.close()
    badkeys.flush()
    badkeys.close()
    return kmers_trie

def analyze_guides(name):
    """Analyze k-mers and find all candidate guideRNAs and their off-targets.

    Load project arguments, build and analyze a trie, find guideRNAs.
    Run after all files were generated by kmers.extract_process_kmers()

    Produce files:
    trie with all k-mers, values store label for good or bad candidate
        guideRNA and coordinates in the genome: <name>/<name>_trie.dat
    intermediate files with candidate guideRNA k-mers used as keys
        to the trie: <name>/<name>_triekeys_v?.txt.gz
    final list of guideRNAs: <name>/<name>_guides.txt.gz

    Args:
    name: project name, used to get project args and in all output

    Return:
    trie.trie object with all k-mers, their coordinates in the genome,
    and labels of good and bad candidate guideRNAs
    """
    # parts = 256

    util.print_log('start analyze_guides()')
    util.print_log('load arguments...')
    args = util.load_args(name)
    util.print_args(args)
    util.print_log('done')
    n = args['greateroffdist']
    parts = 4 ** n

    if os.path.exists('%s%s' % (name,'/blacklist')):
        util.print_log('blacklist directory already exists \n')
        pass
    else:
        os.mkdir(('%s%s' % (name,'/blacklist')))
        util.print_log('blacklist directory made \n')

    # in order to store classified files
    if os.path.exists('%s%s' % (name,'/classifiedfiles')):
        util.print_log('classifiedfiles directory already exists \n')
        pass
    else:
        os.makedirs(('%s%s' % (name,'/classifiedfiles/kmers')))
        os.makedirs(('%s%s' % (name,'/classifiedfiles/triekeys_v1')))
        os.makedirs(('%s%s' % (name,'/classifiedfiles/triekeys_v2')))
        os.makedirs(('%s%s' % (name,'/classifiedfiles/guides')))
        os.makedirs(('%s%s' % (name,'/classifiedfiles/tempfiles')))
        util.print_log('classifiedfiles directory made \n')

    if os.path.exists('%s%s' % (name,'/kmers_tries')):
        util.print_log('kmers_tries directory already exists \n')
        pass
    else:
        os.mkdir(('%s%s' % (name,'/kmers_tries')))
        util.print_log('kmers_tries directory made \n')



    util.print_log('construct trie...')
    kmers_filename = '%s/%s_kmers_shuffled.txt.gz' % (name, name)
    util.print_log('load k-mers from %s' % kmers_filename)
    genome = args['genome']
    goodkeysfile = '%s/%s_triekeys_v1.txt.gz' % (name, name) \
                   if args['altpam'] else ''
    badkeysfile = '%s/%s/%s_nonCandidate_triekeys_with_altpams.txt.gz' % (name,'blacklist',name) if args['altpam'] else ''
    if goodkeysfile:
        util.print_log('print candidate guideRNAs to %s' % goodkeysfile)

    tempdir = '%s%s' % (name,'/classifiedfiles/tempfiles')
    triekeys_v1_dir = '%s%s' % (name,'/classifiedfiles/triekeys_v1')
    triekeys_v1_filenames = ['%s/keys%s.txt.gz' % (triekeys_v1_dir, i) for i in range(parts)]
    kmers_dir = '%s%s' % (name,'/classifiedfiles/kmers')
    kmers_filenames = ['%s/kmers%s.txt.gz' % (kmers_dir, i) for i in range(parts)]

    kmers_trie = build_kmers_trie(kmers_filename, genome, name,
                                  altpam=args['altpam'], pampos=args['pampos'],
                                  maxcount=args['maxoffpos'],goodkeysfile=goodkeysfile,
                                  badkeysfile=badkeysfile,tempdir=tempdir, 
                                  triekeys_v1_filenames=triekeys_v1_filenames, 
                                  kmers_filenames=kmers_filenames, processes=args['threads'], n=n, parts=parts)
    util.print_log('done')


    util.print_log('label as bad guideRNAs multimapping k-mers in trie...')
    # keysinputfile = goodkeysfile if goodkeysfile else kmers_filename
    keysoutputfile = '%s/%s_triekeys_v2.txt.gz' % (name, name)
    nonCandidatekeysoutputfile = '%s/%s/%s_nonCandidate_triekeys_targetSites_with_multiple_perfect_hits.txt.gz' %\
                                 (name,'blacklist',name)
    util.print_log('read keys from %s and write to %s'
                   % (triekeys_v1_dir, keysoutputfile))

    triekeys_v2_dir = '%s%s' % (name,'/classifiedfiles/triekeys_v2')
    triekeys_v2_filenames = ['%s/keys%s.txt.gz' % (triekeys_v2_dir, i) for i in range(parts)]

    filter_keys_trie(tempdir, kmers_trie, triekeys_v1_filenames, triekeys_v2_filenames, keysoutputfile,
                     nonCandidatekeysoutputfile, args['threads'], n, parts)
    util.print_log('done')


    util.print_log('assign correct counts to multimapping k-mers in trie...')
    count_filename = '%s/%s_kmers_counts.txt.gz' % (name, name)
    util.print_log('read counts from %s' % count_filename)
    kmers_trie = label_multimapping(kmers_trie, count_filename, n)
    util.print_log('done')


    util.print_log('label as bad guideRNAs k-mers in trie few mismatches away'
                   ' from other k-mers...')
    sim = args['sim'] - 1
    util.print_log('label as bad k-mers with other k-mer at distance <=%s'
                   % sim)
    keysfile = '%s/%s_triekeys_v2.txt.gz' % (name, name)
    util.print_log('read keys from %s' % keysfile)
    filter_trie_mismatch_similarity(tempdir, name, kmers_trie, args['sim'] - 1, triekeys_v2_filenames, args['threads'], n, parts)
    util.print_log('done')


    util.print_log('produce list of good guideRNAs...')
    keysinputfile = keysfile
    keysoutputfile = '%s/%s_guides.txt.gz' % (name, name)
    nonCandidatekeysoutputfile = '%s/%s/%s_nonCandidate_guides_with_mismatch_neighbors.txt.gz' % (name,'blacklist',name)
    util.print_log('read keys from %s and write to %s'
                   % (keysinputfile, keysoutputfile))

    guides_dir = '%s%s' % (name,'/classifiedfiles/guides')
    guides_filenames = ['%s/%s.txt.gz' % (guides_dir, i) for i in range(parts)]

    filter_keys_trie(tempdir, kmers_trie, triekeys_v2_filenames, guides_filenames, keysoutputfile,
                     nonCandidatekeysoutputfile, args['threads'], n, parts)
    util.print_log('done')

    badkeysfiles = ['%s/badkeys%s.txt.gz' % (tempdir, i) for i in range(parts)]
    for i in range(parts):
        if(os.path.exists(badkeysfiles[i])):
            os.remove(badkeysfiles[i])

    util.print_log('save tries...')
    trie_filename = ['%s/%s/%s_trie%s.dat' % (name, 'kmers_tries', name, i) for i in range(parts)]
    # util.print_log('save in file %s' % trie_filename)
    save_trie(kmers_trie, trie_filename, parts)
    util.print_log('done')

    return kmers_trie
