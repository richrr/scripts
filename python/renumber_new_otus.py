import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re

# Usage: python ~/scripts/python/renumber_new_otus.py -i bisa-signif-taxa-pval-0.01-iv-70-remF8cd80k-table_mc80000.biom.txt -c /home/williamslab/qiime_software/python-2.7.3-release/lib/python2.7/site-packages/picrust/data/16S_13_5_precalculated.tab

def main(args):
    # if it would have been that all the original/standard otu ids for the same taxa had same copy numbers, you could have replaced
    # the 'New*' otu ids (from open reference otu picking) with any original/standard otu id for the same taxa.
    # however, a taxa with different otu ids have different copy numbers, so I cannot replace 'New*' otu ids.
    # For details, refer ~/Data/Ping/qiime_analysis/otus/rosana/removeF8/BlockedIndicatorSpeciesAnalysis/picrust
    parser = argparse.ArgumentParser(description='check which taxa names (if any) have different copy numbers using their original/standard otu ids.')
    parser.add_argument('-i', '--infile')
    #parser.add_argument('-o', '--outfile', default="renumbered-non-new-otu-ids.txt")  # output filename
    parser.add_argument('-c', '--copynumbfile')  # file containing the copy number info for each originial (standard) otu id.
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-n', '--otuidcolumn', default=0, type=int) # 0-index based column used for otu ids
    parser.add_argument('-t', '--taxonnamecolumn', default=-1, type=int) # 0-index based column used for taxonomy names
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    #outfile = args.outfile
    copynumb_file = args.copynumbfile
    delim = args.delimiter
    old_otu_ids_col = args.otuidcolumn
    taxon_col = args.taxonnamecolumn

    d_list=list()
    comment_list=list()
 
    for l in open(infile):
        l = l.strip()
        if '#' not in l:
            d_list.append(l)
        else:
            comment_list.append(l)
            if '#OTU ID' in l or 'taxonomy' in l:
                array = l.split(delim)
                if old_otu_ids_col == array.index('#OTU ID'):
                    pass
                else:
                    print array.index('#OTU ID')
                    sys.exit('Default OTU id column number is incorrect. Use -n \n')
                if array[taxon_col] == 'taxonomy': # using index would fail for taxon, koz default is -1
                    pass
                else:
                    print array.index('taxonomy')
                    sys.exit('Default taxonomy column number is incorrect. Use -t \n')
                    

    taxon_otuid_dict = dict() # dict with taxon as key and (original) otu ids as value
    # print required info for topk
    for l in d_list:
         #print l
         # split contents of line
         array = l.split(delim)
         t = array[taxon_col]
         o = array[old_otu_ids_col]
         if t not in taxon_otuid_dict:
             taxon_otuid_dict[t] = list()
         taxon_otuid_dict[t].append(o)

    otuid_copy_numb_dict = list_to_dict(read_file(copynumb_file)) # dict with (original) otu ids as key and copy number as value

    taxon_copy_numb_dict = dict() # dict with taxon as key and copy numbers as value
    for t in taxon_otuid_dict:
        if t not in taxon_copy_numb_dict:
            taxon_copy_numb_dict[t] = list()
        otus = taxon_otuid_dict[t]
        for o in otus:
            if o in otuid_copy_numb_dict:
                c = otuid_copy_numb_dict[o]
                taxon_copy_numb_dict[t].append(c)
            else: # these are 'New*' otu ids, so absent in the gg 16s count files
                pass

    count = 0
    for t in taxon_otuid_dict:
        v = set(taxon_copy_numb_dict[t])
        if len(v) != 1:
            print t, v
            count += 1
    print '%d taxa, whose otu ids have different copy numbers' %count


#    writeLIST_to_file(list, outfile)

  

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

