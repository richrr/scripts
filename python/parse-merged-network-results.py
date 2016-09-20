import os
import sys
from utils import * 
import re


#python /nfs1/Morgun_Lab/richrr/scripts/python/parse-merged-network-results.py merged_two_expts_pear_ttest_files_corr_merged-parallel-output.csv-consistent_across_expts_output.csv-consistent_across_expts_and_analysis_output.csv ../AnalysToDoList.txt

#python /nfs1/Morgun_Lab/richrr/scripts/python/parse-merged-network-results.py merged_two_expts_pear_ttest_files_comp_merged-parallel-output.csv-consistent_across_expts_output.csv-consistent_across_expts_and_analysis_output.csv ../AnalysToDoList.txt


infile = sys.argv[1]

analysis_file = sys.argv[2]

inlines = read_file(infile)

analys_dict = dict()
analys_numb = 1
for l in read_file(analysis_file):
    l = l.strip()
    if not l:
        continue
    analys_dict[analys_numb] = l.replace('\t', '---||')
    analys_numb += 1

flag_end_header = "FALSE"

n=0
gene_analysis_dict = dict()
found = '' # the anlaysis number obtained from the string

while n < len(inlines)-1:
    l = inlines[n].strip()
    nextl = inlines[n+1].strip()
    #print n+1 # human readable
    if not l:
        flag_end_header = "FALSE"
        n = n+1
        continue
    if (("Analys" in l) and ("AND" in nextl)):
        flag_end_header = "FALSE"
        n = n+2 # skip over the next iteration since next line is going to be "AND"
        continue
    elif (("Analys" in l) and ("AND" not in nextl)) :
        flag_end_header = "TRUE"
        #n = n+1
    if(flag_end_header == "TRUE"):
        #print n
        if nextl: # if not empty
            
            try:
                found = re.search('Analys (\d+) ', l).group(1)
                #print found, analys_dict[int(found)]
            except AttributeError:
                pass
            #make a hash
            for gene in nextl.split(','):
                if gene in gene_analysis_dict:
                    v = gene_analysis_dict[gene]
                    val = "Analysis:" + found + " " + analys_dict[int(found)]
                    v.append(val)
                    gene_analysis_dict[gene] = v
                else:
                    val = "Analysis:" + found + " " + analys_dict[int(found)] 
                    gene_analysis_dict[gene] = [val]
            #print gene_analysis_dict
            #sys.exit()
    n = n+1    
    

out_str = 'FROM: ' + infile + '\n'
for k,v in gene_analysis_dict.items():
    out_str += k + "\n\t" + "\n\t".join(v) + '\n'



outfilename = "./multi-expt-consistent-output-4-lookometer.txt"

if os.path.isfile(outfilename):
   with open(outfilename, "a") as myfile:
      myfile.write(out_str)
else:
   writeTXT_to_file(out_str, outfilename)

    
    
    
