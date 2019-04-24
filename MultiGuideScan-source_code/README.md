Full CRISPR guideRNA analysis pipeline.

Ensure the following dependencies are present

    :::system
        samtools 1.3.1
        easy_install
        coreutils (shuf)
        rename
        python 2.7
        biopython>=1.66
        pysam==0.8.3
        pyfaidx==0.4.7.1
        bx-python==0.7.3
        
        for Rule Set 2 on-target cutting efficiency scores (this must be installed by user: https://pypi.python.org/pypi/scikit-learn/0.16.1)
        sklearn==0.16.1

To install, run 

    :::system
        python setup.py install

This will install binaries `guidescan_processer`, `guidescan_bamdata`, `guidescan_guidequery`, `guidescan_cutting_efficiency_processer`, `guidescan_cutting_efficiency_processer`.

After installation, use

    :::python
        from guidescan import *

in your python session to import all modules of the package, or use

    :::python
        from guidescan import guidequery

to import a particular module and then use functions from the module

    :::python
        guidequery.query_bam()

For local installation, run something like

    :::system
        python setup.py install --user

and then make sure that the local directory with binaries (such as `$HOME/Library/Python/2.7/bin/`) is available in your PATH.


For more info on full pipeline of guideRNA database construction from genomic sequences run

    :::system
        guidescan_processer -h

For more info on guideRNA database construction with custom parameters using previously computed trie run

    :::system
        guidescan_bamdata -h

For more info on accessing precomputed guideRNA database run

    :::system
        guidescan_guidequery -h

For more info on computing cutting efficiency scores for Cas9 20mer gRNAs

	:::system
		guidescan_cutting_efficiency_processer -h

For more info on computing cutting specificity scores for Cas9 20mer gRNAs

	:::system
		guidescan_cutting_specificity_processer -h
