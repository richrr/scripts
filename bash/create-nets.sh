#!/bin/bash

# /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/jackknife-25-samples/merged/gexpress
# usage: bash ~/Morgun_Lab/richrr/scripts/bash/create-nets.sh /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/jackknife-25-samples/merged/gexpress merged_gexp_from_jackknife_from_p0.2_cp0.05_cf_0.1genes_pe_corr_p1_ ~/Morgun_Lab/richrr/Cervical_Cancer/Data/simple-correlation-analysis-only.txt /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/merged/comp/gexpress/p0.2/merged_gexp_tt_comp_p0.2_FolChMedian_merged-parallel-output.csv 1


# bash --version

BPATH=$1 # base path
MRGDFILE=$2 # merged file that has the correlation info
ANALFILE=$3  # analysis file
FCFILE=$4 # merged file that has the comparison info to be used for fold change
FCANALYSNUMB=$5 # analysis number whose fold change number is to be used from the merged file

### make sure to put the right number in "ACCOUNT=?" in the code below
ACCOUNT=1 # put 0 if the analysis file only consists of correlations. currently put 1 since the first analysis is a comparison, the correlations start from the 2nd analysis
### make sure to put the right number in "foldch_?_consistent" in the code below
foldch=0
### make sure to put the right number in "--combPvalCutoff 0.??" in the code below
combPvalCutoff=0.1

echo "Run the different files on SGE"

#http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash

COUNTER=0

declare -a grades
declare -A anums

#http://tldp.org/LDP/abs/html/string-manipulation.html
while IFS=$'\t' read A B C D; do  # this only takes non-empty values. so "high			correlation" only gives values in A and B
    echo "$A $B $C $D"
    let COUNTER=COUNTER+1
    if [ "$A" != "" ] && [ "$B" == "correlation" ]; then
        A=${A// /\\ }
        grades[$COUNTER]=$A
        anums[$A]="$COUNTER"
    fi
done < "$ANALFILE"


#	for grade in "${grades[@]}"; do    # for each grade
#		cmd="$grade ${anums[$grade]}"
#		echo "$cmd" 
#          done




#declare -a fdrs=( 0.2 0.1 0.07 0.05 0.025 0.01 )
#declare -a fdrs=( 0.4 0.3 0.2 )
declare -a fdrs=( 0.4 )

        
        OFILE="cmds-to-create-nets_.txt"
	echo "# run these on sge" > $OFILE

        cmd="cd $BPATH/p1/"
        echo "$cmd" >> $OFILE
        
        frmula="'{print "'$0 > (''"./per_analysis/analysis-" NR "-consistent-corr-dir.txt")}'"'"
        cmd='awk -v RS="" '"$frmula $MRGDFILE""FolChMedian_merged-parallel-output.csv_pval_1_corr_0_foldch_"$foldch"_consistent_results.csv"

	echo "$cmd" >> $OFILE
		
	
	
	############## create networks #################
	
	for fdr in "${fdrs[@]}"; do    # for each fdr cutoff
	
	  for grade in "${grades[@]}"; do    # for each grade
	    	cmd="cd $BPATH/p1/per_analysis/"
		echo "$cmd" >> $OFILE
	
		ndir="$grade/indiv-pval_0.3_comb-pval_$combPvalCutoff""_comb-fdr_$fdr"
		cmd="mkdir -p $ndir"
		echo "$cmd" >> $OFILE

		cmd="cd $ndir/"
		echo "$cmd" >> $OFILE
		
		CORRANALY=$((${anums[$grade]}-ACCOUNT))
	
		cmd="Rscript /nfs3/PHARM/Morgun_Lab/richrr/scripts/R/create-network.R --file ../../../$MRGDFILE""FolChMedian_merged-parallel-output.csv --group $grade --consistent ../../analysis-"$CORRANALY"-consistent-corr-dir.txt --foldchange $FCFILE --indivPvalCutoff 0.3 --combPvalCutoff $combPvalCutoff --analysisfc Analys\ $FCANALYSNUMB\  --analysiscorr Analys\ "${anums[$grade]}"\  --combFDRCutoff $fdr"
		echo "$cmd" >> $OFILE
	
          done
          
        done     

echo "Set execute permissions using 'chmod +x' and run $OFILE"





