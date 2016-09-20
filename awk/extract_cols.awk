#! /bin/csh -f
#
# https://bima.astro.umd.edu/checker/node22.html
# https://www.biostars.org/p/166527/
##                 	 check if called properly
if ($#argv != 3) then
    echo "Usage: $0 mapping_file otu_table output_file"
    echo "The mapping file contains samples to be extracted in first column"
    echo "The Otu table file contains samples as columns"
    echo "For each sample from the mapping file, extract the matching column from the otu table file"
    goto done
endif
##                  save command line args in variables
set mappingfile=$1
set otutable=$2
set outputfile=$3
set headerfile="tmp_header_file.txt"

echo "#OTU ID" > $headerfile
cut -f 1 $mappingfile | grep -v "#" | sort | uniq >> $headerfile
echo "taxonomy" >> $headerfile

awk -F'\t' 'NR==FNR{arr[$1]++;next}{for(i=1; i<=NF; i++) if ($i in arr){a[i]++;}} { for (i in a) printf "%s\t", $i; printf "\n"}' $headerfile $otutable > $outputfile

##           Labels to jump to exit OK (done) or not OK (error)
done:
 exit 0
error:
 exit 1
