#!/bin/bash

echo "Examples of GuideScan usage"
echo "=========================="
echo
echo "Examples show how to run full database generation pipeline"
echo "for ebola genome starting from FASTA files"
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
echo "will result in generation of a folder 'ebola_all/' with"
echo "database for ebola and all intermediate files (including a k-mer trie)"
echo "that can be used for more detailed analysis or for generation"
echo "of the database with a different set of parameters."
echo "Will also generate log files and example output of querying "
echo "the database."
echo 

echo "Make log and file output directories"
command="mkdir logs_ebola"
echo "$" $command
eval $command
echo
command="mkdir guidescan_output_files_ebola"
echo "$" $command
eval $command
echo

echo `date`
echo "Generate guideRNA database from FASTA sequences"
echo "-----------------------------------------------"
command="guidescan_processer -f chromFa/ebola -n ebola_all -d 3 -k 4 -t 8 >logs_ebola/log-guidescan-processer-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Generate more detailed database using intermediates of previous run"
echo "-------------------------------------------------------------------"
command="guidescan_bamdata -n ebola_all --label offdist4 -d 4 -k 5 -t 8  >logs_ebola/log-guidescan-bamdata-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Insert Rule Set 2 cutting efficiency scores"
echo "-------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/ebola_all -n ebola_all_guides.bam -f ./chromFa/ebola/ebola.fa > logs_ebola/log-guidescan-cutting-efficiency-score-insert-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Transfer Rule Set 2 cutting efficiency scores from original database to more detailed database"
echo "--------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/ebola_all -n ebola_all_guides_offdist4.bam -d2 $dir/ebola_all -n2 cutting_efficiency_scores_added_ebola_all_guides.bam >logs_ebola/log-guidescan-cutting-efficiency-score-transfer-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Insert CFD cutting specificity scores"
echo "-------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_specificity_processer -d $dir/ebola_all -n ebola_all_guides.bam -f ./chromFa/ebola/ebola.fa -k $dir/ebola_all > logs_ebola/log-guidescan-cutting-specificity-score-insert-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Insert CFD cutting specificity scores into database with Rule Set 2 cutting efficiency scores"
echo "--------------------------------------------------------"
dir=`pwd`
command="guidescan_cutting_efficiency_processer -d $dir/ebola_all -n cutting_specificity_scores_added_ebola_all_guides.bam -d2 $dir/ebola_all -n2 newdatabase_w_transfered_efficiency_scores.bam >logs_ebola/log-guidescan-cutting-efficiency-score-transfer-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Query the database for specific genomic region (excel output)"
echo "-------------------------------------------------------"
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target within --header -o guidescan_output_files/ --sort offtargets >logs_ebola/log-guidescan-within-query-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target flanks --flankdistance 100 --off --header -o guidescan_output_files/ --sort offtargets >logs_ebola/log-guidescan-flanks-query-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Query the database for specific genomic region and have GuideScan choose candidate sgRNAs (excel output)"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target within --header -o guidescan_output_files/ --select offtargets >logs_ebola/log-guidescan-within-select-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target flanks --flankdistance 100 --header -o guidescan_output_files/ --select offtargets >logs_ebola/log-guidescan-flanks-select-ebola.txt 2>&1" 
echo "$" $command
eval $command
echo "done"
echo

echo "Query the database for specific genomic region (text output)"
echo "-------------------------------------------------------"
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target within --header -o guidescan_output_files/ --sort offtargets --output_format bed >logs_ebola/log-guidescan-within-query-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target flanks --flankdistance 100 --off --header -o guidescan_output_files/ --sort offtargets --output_format bed >logs_ebola/log-guidescan-flanks-query-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Query the database for specific genomic region and have GuideScan choose candidate sgRNAs (text output)"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target within --header -o guidescan_output_files/ --select offtargets --output_format bed >logs_ebola/log-guidescan-within-select-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam -c chrE:400-600 --target flanks --flankdistance 100 --header -o guidescan_output_files/ --select offtargets --output_format bed >logs_ebola/log-guidescan-flanks-select-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo "Query the database for multiple genomic regions in batch and return all sgRNAs"
echo "---------------------------------------------------------"
command="guidescan_guidequery -b ebola_all/ebola_all_guides.bam --batch ./batch/ebola_example_batch_query.txt  --target within --header -o guidescan_output_files/ --sort offtargets >logs_ebola/log-guidescan-within-select-batch-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo

echo `date`
echo "Generate cpf1 guideRNA database from FASTA sequences"
echo "-----------------------------------------------"
command="guidescan_processer -f chromFa/ebola -n ebola_cpf1_all -d 3 -k 4 -t 8 -p TTTN -a XXXX --pampos start >logs_ebola/log-guidescan-cpf1-processer-ebola.txt 2>&1"
echo "$" $command
eval $command
echo "done"
echo `date`
echo
 
