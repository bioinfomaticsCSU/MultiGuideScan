The purpose of the scripts contained within this directory is to allow for the parallel computation of Rule Set 2 on-target cutting efficiency scores (ces) or CFD cutting specificity scores (tss) for a sgRNA database produced by the GuideScan software. Since these database files can be quite large it benefits the user to parallize the ces process. In order to do this the user should do the following:

1.) Split the final database BAM file using BamFileSplitter.py. This script will take the original BAM file and spilt it by lines into a set of smaller BAM files. The original BAM file will not be destroyed in this process. For specific run parameters type python BamFileSplitter.py -h.

2.) Modify ces_parallel_processer.sh $package_dir, $nopos, $full, $fasta, $bam_dir, $split_name, and $PYTHONPATH. Before running the script on your cluster ensure $filename is correct. This variable should not be a path, but rather an explicit name that is derived from your path. Adjust cut -d "/" -f #- to get the filename where the item you are adjusting is #.

3.) Merge the split processed BAM files back into a single database BAM file using BamFileMerger.py. The original split files are not destroyed in this process. For specific run parameters type python BamFileMerger.py -h 
