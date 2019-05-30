Full CRISPR guideRNA analysis pipeline.

To accommodate different CRISPR endonucleases with different specificity requirements, MultiGuideScan allows users to supply the reference genome sequences and define target sequences by setting specific parameters, such as a protospacer adjacent motif (PAM), the PAM positions relative to target sequences, the length of guide RNAs, the Hamming distance *M* (the minimum distance among guide RNAs to ensure that every guide RNA has a unique targeting site in the genome) and the Hamming distance *Q* (the maximum distance considered between the guide RNAs and corresponding off-target sequences). Non-canonical PAM sequences can also be supplied and be helpful for off-target cutting, because they can be recognized and cut by the CRISPR protein with a certain efficiency.

The datasets can be found from [the UCSC genome browser](http://hgdownload.soe.ucsc.edu).

Usage
----

```
Usage: guidescan_processer [options]   or  python processer.py [options]
  Options:
      -h            default: false
  
    * -f            path to fasta file or folder with fasta files (will use all .fa, .fasta, .fa.gz, .fasta.gz files found in the folder)
    * -n            project name, use in all output (will produce a folder with this name containing intermediate and final files in it)
                    default:myguides
    --minchr        minimum chromosome length to consider, chromosomes in input FASTA that are shorter than this will be excluded from 
    				analysis; simple way to exclude scaffolds unassigned to known chromosomes etc.
                    default:10000
      -c            list names of chromosomes (comma separated) that will be used in analysis, or name of file where this list is stored         
                    default: ''
      -l            desired length of guideRNAs (not including PAM)
                    default:20
      -p            PAM sequence
                    default:NGG
      -a            alternative PAM sequences (separate multiple ones by commas), will not be used in primary guideRNAs, but will be considered
      				in off-targets; all PAM sequences should be mutually exclusive and of the same length
                    default:NAG
    --pampos        position of PAM with respect to guideRNA
                    default:end
                    choices:start, end    
      -s            minimum mismatch similarity between guideRNAs; a candidate guideRNA (with primary PAM) should not have alternative occurences
      				in the genome (with any PAM) with less than this many mismatches (not including PAM)
                    default:2  
      -d            maximum distance to consider from guideRNA to its off-target; off-target is an alternative occurrence (with any PAM) of this
                    guideRNA in the genome at edit distance at most this number (including PAM); currently values larger than 4 are infeasible 
                    for large (e.g., mammalian) genomes, and value 3 will take long time to compute; use -1 if do not want any off-target info in 
                    resulting database (can add it later using bamdata)
                    default:3
      -k            a number greater than offdist used for preprocessed data (the length of key for classifying guide RNAs)
                    default:4
    --maxoffpos     maximum number of positions of k-mers to remember; for k-mer occurring multiple times in the genome  (such k-mers cannot be 
                    guideRNAs, but their positions can be off-targets of guideRNAs) store at most this many arbitrary their occurrences in the genome
                    default:10
    --maxoffcount   maximum number of off-targets to store for a guideRNA in a resulting BAM library
                    default:1000
    --gnupath       path to gnu utilities, e.g. "/usr/local/bin"; if empty, use system defaults; requires: cut, sort, uniq, shuf
    * -t            how many processes to use; do not specify more than you have on your system
                    default:1
```

For example:

    python processer.py -f ./chromosomes -n sacCer3_all -d 3 -k 4 -t 10 >logs3/log-guidescan-processer-sacCer3-d3-k4-10.txt 2>&1

    python processer.py -f ./chromosomes -n sacCer3_all -d 4 -k 5 -t 32 >logs3/log-guidescan-processer-sacCer3-d4-k5-32.txt 2>&1
