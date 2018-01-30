import sys
import os
from utils import *
import os.path

#print(sys.argv)

## usage:
# only maps the 0-index
#python /nfs3/PHARM/Morgun_Lab/richrr/scripts/python/map-uniqid-to-symbol-brb.py build_netw_FolChMedian_FCAnalys_12_CorrAnalys_12___li_T2_HFHSNCD_indiv-pval_0.3_comb-pval_0.05_comb-fdr_1_.csv /nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/geneid_symbol_mouse_map.csv 0 1 build_netw_FolChMedian_FCAnalys_12_CorrAnalys_12___li_T2_HFHSNCD_indiv-pval_0.3_comb-pval_0.05_comb-fdr_1_.csv-gene-symb.csv network

# maps first two columns (See last arg)
#python /nfs3/PHARM/Morgun_Lab/richrr/scripts/python/map-uniqid-to-symbol-brb.py build_netw_FolChMedian_FCAnalys_12_CorrAnalys_12___li_T2_HFHSNCD_indiv-pval_0.3_comb-pval_0.05_comb-fdr_1_.csv /nfs3/PHARM/Morgun_Lab/richrr/db/TopHat/geneid_symbol_mouse_map.csv 0 1 build_netw_FolChMedian_FCAnalys_12_CorrAnalys_12___li_T2_HFHSNCD_indiv-pval_0.3_comb-pval_0.05_comb-fdr_1_.csv-gene-symb.csv network 0,1



file1 = sys.argv[1]
file2 = sys.argv[2]
colf1 = int(sys.argv[3])
colf2 = int(sys.argv[4])
ofile = sys.argv[5]

delim = ","
joiner = ","
col_to_change = '0'
replace_suffix = ''

# first check if the new map id col behaves correctly for 0 index change, similar to its predecessor function
# then add more columns as arguments to this code and see if it works

# to move the first partner1<==>partner2 column to end
def sanitize_net_file(file1):
  line1 = read_file(file1)[0] # treated as a list
  line1 = line1.strip()
  #print line1
  conts = line1.split(delim)
  #print conts
  
  # output dict needs a list for new column ordering
  # move the first column id to last
  fieldnames = conts[1:] + conts[:1]
  #print fieldnames
  
  #https://stackoverflow.com/questions/33001490/python-re-ordering-columns-in-a-csv
  import csv
  sanfile = 'reordered-' + file1
  if os.path.exists(sanfile):
      sys.exit("reordered file exists! Delete before proceeding")
  with open(file1, 'r') as infile, open(sanfile, 'a') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    # reorder the header first
    writer.writeheader()
    for row in csv.DictReader(infile):
        # writes the reordered rows to the new file
        writer.writerow(row)
    return sanfile



if(len(sys.argv) > 6): # since the sys.argv[0] is the script name.
    #delim = sys.argv[6]
    #joiner = sys.argv[7]
    if(sys.argv[6] == "network"):
        file1 = sanitize_net_file(file1)


if(len(sys.argv) > 7): # since the sys.argv[0] is the script name.
    col_to_change = sys.argv[7]


if(len(sys.argv) > 8):
    replace_suffix = sys.argv[8]



# only maps the 0-index
#ret_list = map_identifiers(file1, file2, delim , joiner, colf1, colf2) 

# can map more than one column
ret_list = map_identifiers_col(file1, file2, col_to_change, delim , joiner, colf1, colf2, missingidasNA=False, replace_suffix = replace_suffix) 
writeLIST_to_file(ret_list, ofile)

