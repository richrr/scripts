#!/usr/bin/env gawk

###
#Robert Settlage, PhD
#Data Analysis Core
#August 2013
###

###needs pass in "sample_name read_prefix" "header as per qiime format "orig_bc=GGAGACAAGGGA new_bc=GGAGACAAGGGA bc_diffs=0" and file to convert
###converts to a fasta
BEGIN {
        sample = ARGV[1]
        header = ARGV[2]
        i=0
        j=0
        while (getline < ARGV[3]) {
	   i++
	   if (i==1) {
              print sample"_"j" "header
	      j++
	   }else if(i==2){
	      print $0
           }
	   if(i==4){
	      i=0
	   q}
        }
        exit (0)
}
