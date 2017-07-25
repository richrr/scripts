#! /bin/bash

#https://gist.github.com/slavailn/a139c2414983334ea2d33daade49db7b/


indir=$1

echo -e "SampleID\ttotal_reads\tmapped_reads\tmultiple_alignments\tmapped%\tunique_reads\tunique_mapped%"

for dir in $indir/*
do

    dirname=`basename $dir`
    samplename=${dirname%_tophat_out}
    cd $dir
    
    total=""
    mapped=""
    multiple=""
    perc_mapped=""
    unique=""
    perc_unique_mapped=""
    
    while read line 
    do 
         if [[ $line == *"Input"* ]]; then
             total=`echo $line | awk '{print $3}'`
             
         fi

         if [[ $line == *"Mapped"* ]]; then
             mapped=`echo $line | awk '{print $3}'`
             
         fi

         if [[ $line == *"of these"* ]]; then
             multiple=`echo $line | awk '{print $3}'`
             
         fi          

    done < align_summary.txt
    
    perc_mapped=`echo $mapped*100/$total | bc`
    unique=`echo $mapped-$multiple | bc`
    perc_unique_mapped=`echo $unique*100/$total | bc`
    echo -e "$samplename\t$total\t$mapped\t$multiple\t$perc_mapped\t$unique\t$perc_unique_mapped"

done
