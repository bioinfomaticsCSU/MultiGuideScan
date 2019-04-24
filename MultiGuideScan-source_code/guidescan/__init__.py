"""Full CRISPR-Cas guideRNA analysis pipeline.

Submodules:
processer: Master script for guideRNA database construction
kmers: Extract k-mers from genomic sequence and do initial processing
guides: Analyze k-mers and extract guideRNAs with all relevant info about them
bamdata: Create and manipulate guideRNA database in BAM format
guidequery: User interface script to precomputed guideRNA database
util: Utilities
trie: trie implementation from Biopython, modified

Before use, run the following commands to compile trie module:
$ cd trie
$ python setup.py build
$ cd ..
$ ln -s trie/build/lib*/trie.so

For more info on full pipeline of guideRNA database construction
from genomic sequences run

$ python processer.py -h

For more info on guideRNA database construction with custom parameters using
previously computed trie run

$ python bamdata.py -h

For more info on accessing precomputed guideRNA database run

$ python guidequery.py -h
"""

__all__ = ['util', 'kmers', 'guides', 'bamdata', 'guidequery','ces_compute','ces_insert','ces_database_to_database_transfer']
