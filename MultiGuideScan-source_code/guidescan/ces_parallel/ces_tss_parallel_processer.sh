#!/bin/bash

set -e

package_dir="absolute-filepath-to-GuideScan-package"
nopos="absolute-filepath-to-/V3_model_nopos.pickle"
full="absolute-filepath-to-/V3_model_full.pickle"
fasta="absolute-filepath-to-organism-fasta-file-and-its-index/ebola.fa"
bam_dir="absolute-filepath-to-directory-hosting-split-bam-files"
split_name="split-file-prefix-name"

counter=1

export PYTHONPATH="python-PATH":$PYTHONPATH

for FILE in $(find $bam_dir -name $split_name*)

	do 
		output1="ces_tss_parallel_"$filename
		filename=$(echo $FILE | cut -d "/" -f 11-)
		#echo $filename
		echo "------------------------------"
		echo "	ces/tss parallel " $filename
		echo "------------------------------"
		echo 
		#echo "---cutting efficiency---"
		#echo
		#echo "python $package_dir/cutting_efficiency_processer.py -d $bam_dir -n $filename -f $fasta --nopos $nopos --full $full" #| qsub -l nodes=1:ppn=1,mem=1gb,walltime=00:00:20 -N $output1
		#echo "---cutting specificity---"
		#echo
		#echo "python $package_dir/cutting_specificity_processer.py -d $bam_dir -n $filename -f $fasta" #| qsub -l nodes=1:ppn=1,mem=1gb,walltime=00:00:20 -N $output1
		counter=$((counter+1))
	done  

