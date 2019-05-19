#!/bin/bash

echo "Examples of GuideScan usage"
echo "=========================="
echo
echo "Examples show how to run full database generation pipeline"
echo "for yeast genome starting from FASTA files"
echo "and how to query the resulting database to get info"
echo "about guideRNAs in specific genomic regions."
echo
echo "Assumes uniguide package was already installed."
echo "FASTA files were downloaded from"
echo "http://hgdownload.cse.ucsc.edu/"
echo
echo "This is script guidescan-example.sh"
echo "Uncomment any lines in the script starting with 'eval' to run"
echo "specific examples. Uncommenting all examples and using them all"
echo "will result in generation of a folder 'yeast_all/' with"
echo "database for yeast and all intermediate files (including a k-mer trie)"
echo "that can be used for more detailed analysis or for generation"
echo "of the database with a different set of parameters."
echo "Will also generate log files and example output of querying "
echo "the database."
echo 

echo `date`
echo "Make log and file output directories"
command="mkdir logs"
echo "$" $command
eval $command
echo
command="mkdir guidescan_output_files"
echo "$" $command
eval $command
echo

echo "Generate guideRNA database from FASTA sequences"
echo "-----------------------------------------------"
command="guidescan_processer -f chromFa/chromosomes -n yeast_all -d 3 -k 4 -t 10 >logs/log-guidescan-processer-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Generate more detailed database using intermediates of previous run"
echo "-------------------------------------------------------------------"
command="guidescan_bamdata -n yeast_all --label offdist4 -d 4 -k 5 -t 10 >logs/log-guidescan-bamdata-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Insert Rule Set 2 cutting efficiency scores"
echo "-------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/yeast_all -n yeast_all_guides.bam -f ./chromFa/cutting_efficiency/yeast.fa > logs/log-guidescan-cutting-efficiency-score-insert-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Transfer Rule Set 2 cutting efficiency scores from original database to more detailed database"
echo "--------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/yeast_all -n yeast_all_guides_offdist4.bam -d2 $dir/yeast_all -n2 cutting_efficiency_scores_added_yeast_all_guides.bam >logs/log-guidescan-cutting-efficiency-score-transfer-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Insert CFD cutting specificity scores"
echo "-------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_specificity_processer -d $dir/yeast_all -n yeast_all_guides.bam -f ./chromFa/cutting_efficiency/yeast.fa -k $dir/yeast_all > logs/log-guidescan-cutting-specificity-score-insert-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Insert CFD cutting specificity scores into database with Rule Set 2 cutting efficiency scores"
echo "--------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/yeast_all -n cutting_specificity_scores_added_yeast_all_guides.bam -d2 $dir/yeast_all -n2 newdatabase_w_transfered_efficiency_scores.bam >logs/log-guidescan-cutting-specificity-efficiency-score-transfer-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Query the database for specific genomic region (excel output)"
echo "-------------------------------------------------------"
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target within --header -o guidescan_output_files/ --sort offtargets >logs/log-guidescan-within-query-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target flanks --flankdistance 100 --off --header -o guidescan_output_files/ --sort offtargets >logs/log-guidescan-flanks-query-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Query the database for specific genomic region and have GuideScan choose candidate sgRNAs (excel output)"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target within --header -o guidescan_output_files/ --select offtargets >logs/log-guidescan-within-select-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target flanks --flankdistance 100 --header -o guidescan_output_files/ --select offtargets >logs/log-guidescan-flanks-select-yeast.txt 2>&1" 
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Query the database for specific genomic region (text output)"
echo "-------------------------------------------------------"
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target within --header -o guidescan_output_files/ --sort offtargets --output_format bed >logs/log-guidescan-within-query-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target flanks --flankdistance 100 --off --header -o guidescan_output_files/ --sort offtargets --output_format bed >logs/log-guidescan-flanks-query-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Query the database for specific genomic region and have GuideScan choose candidate sgRNAs (text output)"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target within --header -o guidescan_output_files/ --select offtargets --output_format bed >logs/log-guidescan-within-select-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam -c chrI:4000-6000 --target flanks --flankdistance 100 --header -o guidescan_output_files/ --select offtargets --output_format bed >logs/log-guidescan-flanks-select-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Query the database for multiple genomic regions in batch and return all sgRNAs"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b yeast_all/yeast_all_guides.bam --batch ./batch/yeast_example_batch_query.txt  --target within --header -o guidescan_output_files/ --sort offtargets >logs/log-guidescan-within-select-batch-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Generate cpf1 guideRNA database from FASTA sequences"
echo "-----------------------------------------------"
command="guidescan_processer -f chromFa/chromosomes -n yeast_cpf1_all -d 3 -k 4 -t 10 -p TTTN -a XXXX --pampos start >logs/log-guidescan-cpf1-processer-yeast.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo `date`
echo
 
