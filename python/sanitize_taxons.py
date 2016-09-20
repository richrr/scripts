# to sanitize the taxonomy names from QIIME
   
import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re

# Usage: python ~/scripts/python/sanitize_taxons.py -i summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt -o sanitized_summary_otu_table_mc2_w_tax_no_pynast_failures.biom.txt-rep_set_tax_assign.txt
def main(args):

    parser = argparse.ArgumentParser(description='Sanitize taxonomy names/OTU.')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="out.txt")  # output filename
    parser.add_argument('-s', '--sanitizecolumn', type =int, default=1)  # 0-index based column to be sanitized
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfile = args.outfile
    delim = args.delimiter
    sanitcol = args.sanitizecolumn


    # k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Kaistobacter; s__
    taxonomy_index = {'k' : 0, 'p' : 1, 'c' : 2, 'o' : 3, 'f' : 4, 'g' : 5, 's' : 6}
    # list of tuples, ordered as per the index/value of the dictionary
    sorted_taxonomy_index = sorted(taxonomy_index.iteritems(), key=operator.itemgetter(1))


    lines = read_file(infile)
    outlist = list()
    for l in lines:
        l = l.strip()
        if '#' in l:
            outlist.append(l)
        else:
            cont = l.split(delim)
            sanitized_value = cont[sanitcol]
            if sanitized_value == 'Unassigned':
                sanitized_value = 'k__Unassigned; p__; c__; o__; f__; g__; s__'
            elif sanitized_value == 'None':
                sanitized_value = 'k__None; p__; c__; o__; f__; g__; s__'
            elif sanitized_value == 'Unclassified':
                sanitized_value = 'k__Unclassified; p__; c__; o__; f__; g__; s__'
            else:
                for tup in sorted_taxonomy_index:
                    taxon=tup[0]+'__'
                    if taxon not in sanitized_value:
                        sanitized_value = sanitized_value + '; ' + taxon
                sanitized_value = sanitized_value.strip('; ')
            cont[sanitcol] = sanitized_value
            outlist.append(delim.join(cont))
    writeLIST_to_file(outlist, outfile)

  

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

